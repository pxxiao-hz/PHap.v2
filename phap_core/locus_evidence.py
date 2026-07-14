"""Evaluate complete PAF input into auditable mT2T locus-routing decisions."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Optional, Tuple

from .alignment_chains import PafChainEvidence, summarize_paf_chains
from .locus_rescue import (
    LocusRescueDecision,
    LocusRescueThresholds,
    decide_locus_rescue,
)
from .paf import PafRecord, PafRecordDecision, select_primary_records


@dataclass(frozen=True)
class PafAlignmentAudit:
    """Final filtering status of one parsed PAF record."""

    source: str
    line_number: int
    unitig_id: str
    locus_id: str
    alignment_type: Optional[str]
    status: str
    reason: str
    final_status: str
    final_reason: str


@dataclass(frozen=True)
class LocusEvidenceResult:
    """All deterministic products of strict PAF locus evaluation."""

    allowed_loci: Tuple[str, ...]
    alignment_audit: Tuple[PafAlignmentAudit, ...]
    chains: Tuple[PafChainEvidence, ...]
    decisions: Tuple[LocusRescueDecision, ...]
    filtered_records: Tuple[PafRecord, ...]


def evaluate_locus_evidence(
    records: Sequence[PafRecord],
    *,
    unitig_ids: Iterable[str],
    source_states: Mapping[str, str],
    independent_read_support: Mapping[str, int],
    allowed_loci: Iterable[str],
    thresholds: LocusRescueThresholds,
    min_query_length: int,
    min_alignment_block_length: int,
    max_query_gap: int,
    max_target_gap: int,
    require_tp: bool = True,
) -> LocusEvidenceResult:
    """Evaluate full PAF evidence before any best-locus information is lost."""

    if min_query_length < 1:
        raise ValueError("min_query_length must be positive")
    if min_alignment_block_length < 1:
        raise ValueError("min_alignment_block_length must be positive")
    loci = tuple(sorted(set(allowed_loci)))
    if any(not locus for locus in loci):
        raise ValueError("allowed locus IDs must not be empty")
    locus_set = set(loci)
    all_unitigs = set(unitig_ids)
    if any(not unitig_id for unitig_id in all_unitigs):
        raise ValueError("unitig IDs must not be empty")
    all_unitigs.update(record.query_name for record in records)
    for unitig_id, count in independent_read_support.items():
        if not unitig_id or isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ValueError(
                "independent read support must map non-empty unitig IDs "
                "to non-negative integers"
            )
        all_unitigs.add(unitig_id)

    _, primary_audit = select_primary_records(
        records,
        require_tp=require_tp,
    )
    primary_by_line = {
        (row.source, row.line_number): row for row in primary_audit
    }
    accepted_for_chains = []
    alignment_audit = []
    seen_alignment_signatures: set[Tuple[object, ...]] = set()
    for record in sorted(records, key=lambda row: (row.source, row.line_number)):
        primary = primary_by_line[(record.source, record.line_number)]
        status, reason = _record_filter_status(
            record,
            primary,
            locus_set=locus_set,
            min_query_length=min_query_length,
            min_alignment_block_length=min_alignment_block_length,
        )
        if status == "accepted":
            signature = _alignment_signature(record)
            if signature in seen_alignment_signatures:
                status, reason = "excluded", "duplicate_alignment"
            else:
                seen_alignment_signatures.add(signature)
                accepted_for_chains.append(record)
        alignment_audit.append(
            PafAlignmentAudit(
                source=record.source,
                line_number=record.line_number,
                unitig_id=record.query_name,
                locus_id=record.target_name,
                alignment_type=record.alignment_type,
                status=status,
                reason=reason,
                final_status="pending" if status == "accepted" else "excluded",
                final_reason="pending_locus_decision" if status == "accepted" else reason,
            )
        )

    chains = summarize_paf_chains(
        accepted_for_chains,
        max_query_gap=max_query_gap,
        max_target_gap=max_target_gap,
    )
    chains_by_unitig = defaultdict(list)
    for chain in chains:
        chains_by_unitig[chain.unitig_id].append(chain)

    decisions = []
    for unitig_id in sorted(all_unitigs):
        unitig_chains = tuple(
            sorted(
                chains_by_unitig.get(unitig_id, ()),
                key=lambda chain: chain.locus_id,
            )
        )
        decisions.append(
            decide_locus_rescue(
                unitig_id,
                [chain.as_locus_summary() for chain in unitig_chains],
                source_state=source_states.get(unitig_id, "missing_dosage"),
                independent_read_support=independent_read_support.get(unitig_id, 0),
                thresholds=thresholds,
            )
        )

    assigned_loci = {
        decision.unitig_id: decision.assigned_locus
        for decision in decisions
        if decision.assigned_locus is not None
    }
    accepted_record_keys = {
        (record.source, record.line_number)
        for chain in chains
        if assigned_loci.get(chain.unitig_id) == chain.locus_id
        for record in chain.records
    }
    assigned_loci_or_none = {
        decision.unitig_id: decision.assigned_locus for decision in decisions
    }
    final_alignment_audit = []
    for row in alignment_audit:
        if row.status != "accepted":
            final_alignment_audit.append(row)
            continue
        record_key = (row.source, row.line_number)
        assigned_locus = assigned_loci_or_none.get(row.unitig_id)
        if record_key in accepted_record_keys:
            final_status = "written"
            final_reason = "selected_assigned_locus"
        elif assigned_locus is None:
            final_status = "excluded"
            final_reason = "unitig_unassigned"
        else:
            final_status = "excluded"
            final_reason = "competing_locus_not_selected"
        final_alignment_audit.append(
            PafAlignmentAudit(
                source=row.source,
                line_number=row.line_number,
                unitig_id=row.unitig_id,
                locus_id=row.locus_id,
                alignment_type=row.alignment_type,
                status=row.status,
                reason=row.reason,
                final_status=final_status,
                final_reason=final_reason,
            )
        )
    filtered_records = tuple(
        sorted(
            (
                record
                for record in accepted_for_chains
                if (record.source, record.line_number) in accepted_record_keys
            ),
            key=lambda record: (
                record.query_name,
                record.target_name,
                record.query_start,
                record.query_end,
                record.target_start,
                record.target_end,
                record.line_number,
            ),
        )
    )
    return LocusEvidenceResult(
        allowed_loci=loci,
        alignment_audit=tuple(final_alignment_audit),
        chains=chains,
        decisions=tuple(decisions),
        filtered_records=filtered_records,
    )


def _record_filter_status(
    record: PafRecord,
    primary: PafRecordDecision,
    *,
    locus_set: set[str],
    min_query_length: int,
    min_alignment_block_length: int,
) -> Tuple[str, str]:
    if primary.status != "accepted":
        return primary.status, primary.reason
    if record.query_length < min_query_length:
        return "excluded", "query_below_min_length"
    if record.alignment_block_length < min_alignment_block_length:
        return "excluded", "alignment_below_min_block_length"
    if record.target_name not in locus_set:
        return "excluded", "target_not_selected"
    return "accepted", "chain_candidate"


def _alignment_signature(record: PafRecord) -> Tuple[object, ...]:
    return (
        record.query_name,
        record.query_length,
        record.query_start,
        record.query_end,
        record.strand,
        record.target_name,
        record.target_length,
        record.target_start,
        record.target_end,
        record.matching_bases,
        record.alignment_block_length,
        record.mapping_quality,
    )
