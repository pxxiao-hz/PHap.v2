"""Conflict-aware oriented overlap graphs for mT2T sequence construction."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Dict, Optional, Tuple

from .alignment_chains import PafChainEvidence, summarize_paf_chains
from .fasta import FastaRecord
from .locus_rescue import interval_union_length
from .mt2t_graph import (
    _extract_paths,
    _final_alignment_reasons,
    _materialize_outputs,
    _resolve_graph,
)
from .mt2t_model import (
    _ClassifiedChain,
    AlignmentAudit,
    ChainAudit,
    ContainmentAudit,
    EdgeAudit,
    Mt2tAssemblyResult,
    Mt2tOverlapParameters,
    OrientedContig,
    OrientedOverlapEdge,
    PairChainEvidence,
    reverse_complement,
)
from .paf import PafRecord, select_primary_records

__all__ = [
    "Mt2tAssemblyResult",
    "Mt2tOverlapParameters",
    "assemble_mt2t",
    "reverse_complement",
]


def assemble_mt2t(
    fasta_records: Sequence[FastaRecord],
    paf_records: Sequence[PafRecord],
    parameters: Mt2tOverlapParameters,
) -> Mt2tAssemblyResult:
    """Build mT2T paths while preserving every input contig in routing."""

    records_by_id = {record.identifier: record for record in fasta_records}
    if len(records_by_id) != len(fasta_records):
        raise ValueError("FASTA identifiers must be unique")
    sequence_lengths = {
        identifier: len(record.sequence) for identifier, record in records_by_id.items()
    }
    _validate_paf_fasta_boundary(paf_records, sequence_lengths)
    candidate_records, initial_audit = _filter_paf_records(
        paf_records,
        parameters,
    )
    raw_chains = summarize_paf_chains(
        candidate_records,
        max_query_gap=parameters.max_query_gap,
        max_target_gap=parameters.max_target_gap,
    )
    pair_chains = tuple(_pair_chain(chain) for chain in raw_chains)
    classified = tuple(
        _classify_chain(chain, parameters) for chain in pair_chains
    )

    edge_candidates = tuple(
        item.edge for item in classified if item.edge is not None
    )
    selected_edges, duplicate_edge_audits = _deduplicate_edges(edge_candidates)
    contained, containment_audit = _resolve_containment(
        classified,
        selected_edges,
    )
    graph_edges = tuple(
        edge
        for edge in selected_edges
        if edge.left.contig_id not in contained
        and edge.right.contig_id not in contained
    )
    accepted_edges, graph_edge_audits, node_conflicts = _resolve_graph(graph_edges)
    for row in containment_audit:
        if row.reason in {
            "competing_containment_hosts",
            "terminal_overlap_conflicts_with_containment",
        }:
            node_conflicts.setdefault(row.child_id, row.reason)
    paths = _extract_paths(accepted_edges)

    outputs, routing, joins = _materialize_outputs(
        fasta_records,
        paths,
        contained,
        node_conflicts,
        parameters,
    )
    final_reason_by_line = _final_alignment_reasons(
        classified,
        duplicate_edge_audits + graph_edge_audits,
        containment_audit,
    )
    alignment_audit = tuple(
        AlignmentAudit(
            source=row.source,
            line_number=row.line_number,
            query_id=row.query_id,
            target_id=row.target_id,
            initial_status=row.initial_status,
            initial_reason=row.initial_reason,
            final_status=(
                "evaluated" if row.initial_status == "accepted" else "excluded"
            ),
            final_reason=(
                final_reason_by_line.get((row.source, row.line_number), "not_selected")
                if row.initial_status == "accepted"
                else row.initial_reason
            ),
        )
        for row in initial_audit
    )
    chain_audit = tuple(
        ChainAudit(
            query_id=item.chain.query_id,
            target_id=item.chain.target_id,
            strand=item.chain.strand,
            record_count=item.chain.record_count,
            identity=item.chain.identity,
            query_coverage=item.chain.query_coverage,
            target_coverage=item.chain.target_coverage,
            query_union_bases=item.chain.query_union_bases,
            target_union_bases=item.chain.target_union_bases,
            query_span_bases=item.chain.query_end - item.chain.query_start,
            target_span_bases=item.chain.target_end - item.chain.target_start,
            score=item.chain.score,
            collinear=item.chain.collinear,
            classification=item.classification,
            reason=item.reason,
        )
        for item in sorted(
            classified,
            key=lambda row: (row.chain.query_id, row.chain.target_id),
        )
    )
    return Mt2tAssemblyResult(
        alignment_audit=alignment_audit,
        chain_audit=chain_audit,
        edge_audit=tuple(
            sorted(
                duplicate_edge_audits + graph_edge_audits,
                key=lambda row: row.edge.key,
            )
        ),
        containment_audit=containment_audit,
        paths=paths,
        outputs=outputs,
        routing=routing,
        joins=joins,
        contained_ids=tuple(sorted(contained)),
    )



def _filter_paf_records(
    records: Sequence[PafRecord],
    parameters: Mt2tOverlapParameters,
) -> Tuple[Tuple[PafRecord, ...], Tuple[AlignmentAudit, ...]]:
    _, primary_audit = select_primary_records(records, require_tp=True)
    primary_by_key = {
        (row.source, row.line_number): row for row in primary_audit
    }
    accepted: list[PafRecord] = []
    audits: list[AlignmentAudit] = []
    signatures = set()
    for record in sorted(records, key=lambda row: (row.source, row.line_number)):
        primary = primary_by_key[(record.source, record.line_number)]
        status = primary.status
        reason = primary.reason
        if status == "accepted" and record.query_name == record.target_name:
            status, reason = "excluded", "self_alignment"
        elif (
            status == "accepted"
            and record.alignment_block_length
            < parameters.min_alignment_block_length
        ):
            status, reason = "excluded", "alignment_below_min_block_length"
        elif status == "accepted" and (
            record.query_length < parameters.min_contig_length
            or record.target_length < parameters.min_contig_length
        ):
            status, reason = "excluded", "contig_below_graph_min_length"
        if status == "accepted":
            signature = _record_signature(record)
            if signature in signatures:
                status, reason = "excluded", "duplicate_alignment"
            else:
                signatures.add(signature)
                accepted.append(record)
                reason = "chain_candidate"
        audits.append(
            AlignmentAudit(
                source=record.source,
                line_number=record.line_number,
                query_id=record.query_name,
                target_id=record.target_name,
                initial_status=status,
                initial_reason=reason,
                final_status="pending" if status == "accepted" else "excluded",
                final_reason="pending_chain_evaluation" if status == "accepted" else reason,
            )
        )
    return tuple(accepted), tuple(audits)


def _pair_chain(chain: PafChainEvidence) -> PairChainEvidence:
    target_union = interval_union_length(
        chain.target_intervals,
        sequence_length=chain.target_length,
    )
    target_coverage = target_union / chain.target_length
    return PairChainEvidence(
        query_id=chain.unitig_id,
        target_id=chain.locus_id,
        strand=chain.strand,
        query_length=chain.query_length,
        target_length=chain.target_length,
        query_start=min(start for start, _ in chain.query_intervals),
        query_end=max(end for _, end in chain.query_intervals),
        target_start=min(start for start, _ in chain.target_intervals),
        target_end=max(end for _, end in chain.target_intervals),
        query_coverage=chain.query_coverage,
        target_coverage=target_coverage,
        query_union_bases=chain.union_query_bases,
        target_union_bases=target_union,
        identity=chain.identity,
        score=chain.identity * min(chain.query_coverage, target_coverage),
        collinear=chain.collinear,
        record_count=chain.record_count,
        records=chain.records,
    )


def _classify_chain(
    chain: PairChainEvidence,
    parameters: Mt2tOverlapParameters,
) -> _ClassifiedChain:
    if not chain.collinear:
        return _ClassifiedChain(chain, "rejected", "non_collinear_chain")
    if chain.identity < parameters.min_identity:
        return _ClassifiedChain(chain, "rejected", "identity_below_threshold")
    shorter_coverage, longer_coverage = _short_long_coverages(chain)
    if shorter_coverage < parameters.min_shorter_coverage:
        return _ClassifiedChain(chain, "rejected", "shorter_coverage_below_threshold")
    if longer_coverage < parameters.min_longer_coverage:
        return _ClassifiedChain(chain, "rejected", "longer_coverage_below_threshold")

    containment = _containment_candidate(chain, parameters)
    if containment is not None:
        child_id, host_id = containment
        return _ClassifiedChain(
            chain,
            "containment_candidate",
            "short_contig_covered_inside_long_contig",
            child_id=child_id,
            host_id=host_id,
        )

    edge = _terminal_edge(chain)
    if edge is None:
        return _ClassifiedChain(chain, "rejected", "not_a_unique_terminal_overlap")
    return _ClassifiedChain(
        chain,
        "terminal_overlap_candidate",
        "unique_suffix_prefix_geometry",
        edge=edge,
    )


def _containment_candidate(
    chain: PairChainEvidence,
    parameters: Mt2tOverlapParameters,
) -> Optional[Tuple[str, str]]:
    if chain.query_length == chain.target_length:
        return None
    if chain.query_length < chain.target_length:
        child_id, host_id = chain.query_id, chain.target_id
        child_coverage = chain.query_coverage
        host_length = chain.target_length
        host_start, host_end = chain.target_start, chain.target_end
        child_start, child_end = chain.query_start, chain.query_end
        child_length = chain.query_length
    else:
        child_id, host_id = chain.target_id, chain.query_id
        child_coverage = chain.target_coverage
        host_length = chain.query_length
        host_start, host_end = chain.query_start, chain.query_end
        child_start, child_end = chain.target_start, chain.target_end
        child_length = chain.target_length
    margin = int(host_length * parameters.internal_margin_ratio)
    internal = host_start >= margin and host_length - host_end >= margin
    completely_covered = (
        child_coverage == 1.0
        and child_start == 0
        and child_end == child_length
    )
    if completely_covered and internal:
        return child_id, host_id
    return None


def _terminal_edge(
    chain: PairChainEvidence,
) -> Optional[OrientedOverlapEdge]:
    if (
        chain.query_union_bases != chain.query_end - chain.query_start
        or chain.target_union_bases != chain.target_end - chain.target_start
    ):
        return None
    query_orientation = "+" if chain.strand == "+" else "-"
    query_interval = _orient_interval(
        (chain.query_start, chain.query_end),
        chain.query_length,
        query_orientation,
    )
    target_interval = (chain.target_start, chain.target_end)
    candidates = []
    if (
        target_interval[1] == chain.target_length
        and query_interval[0] == 0
    ):
        candidates.append(
            _canonical_edge(
                OrientedContig(chain.target_id, "+"),
                OrientedContig(chain.query_id, query_orientation),
                chain.target_length,
                chain.query_length,
                target_interval,
                query_interval,
                chain,
            )
        )
    if (
        query_interval[1] == chain.query_length
        and target_interval[0] == 0
    ):
        candidates.append(
            _canonical_edge(
                OrientedContig(chain.query_id, query_orientation),
                OrientedContig(chain.target_id, "+"),
                chain.query_length,
                chain.target_length,
                query_interval,
                target_interval,
                chain,
            )
        )
    unique = {candidate.key: candidate for candidate in candidates}
    if len(unique) != 1:
        return None
    return next(iter(unique.values()))


def _canonical_edge(
    left: OrientedContig,
    right: OrientedContig,
    left_length: int,
    right_length: int,
    left_interval: Tuple[int, int],
    right_interval: Tuple[int, int],
    chain: PairChainEvidence,
) -> OrientedOverlapEdge:
    edge = OrientedOverlapEdge(
        left=left,
        right=right,
        left_length=left_length,
        right_length=right_length,
        left_interval=left_interval,
        right_interval=right_interval,
        query_id=chain.query_id,
        target_id=chain.target_id,
        identity=chain.identity,
        query_coverage=chain.query_coverage,
        target_coverage=chain.target_coverage,
        score=chain.score,
        source_records=tuple(
            sorted((record.source, record.line_number) for record in chain.records)
        ),
    )
    reversed_edge = edge.reversed()
    forward_signature = (
        edge.left.contig_id,
        edge.left.orientation,
        edge.right.contig_id,
        edge.right.orientation,
    )
    reverse_signature = (
        reversed_edge.left.contig_id,
        reversed_edge.left.orientation,
        reversed_edge.right.contig_id,
        reversed_edge.right.orientation,
    )
    return reversed_edge if reverse_signature < forward_signature else edge


def _deduplicate_edges(
    edges: Sequence[OrientedOverlapEdge],
) -> Tuple[Tuple[OrientedOverlapEdge, ...], Tuple[EdgeAudit, ...]]:
    grouped: Dict[Tuple[object, ...], list[OrientedOverlapEdge]] = defaultdict(list)
    for edge in edges:
        grouped[edge.key].append(edge)
    selected: list[OrientedOverlapEdge] = []
    audits: list[EdgeAudit] = []
    for key in sorted(grouped):
        ranked = sorted(grouped[key], key=_edge_rank)
        selected.append(ranked[0])
        for duplicate in ranked[1:]:
            audits.append(EdgeAudit(duplicate, "rejected", "duplicate_edge_lower_score"))
    return tuple(selected), tuple(audits)


def _resolve_containment(
    classified: Sequence[_ClassifiedChain],
    edges: Sequence[OrientedOverlapEdge],
) -> Tuple[Dict[str, str], Tuple[ContainmentAudit, ...]]:
    by_child: Dict[str, list[_ClassifiedChain]] = defaultdict(list)
    for item in classified:
        if item.child_id is not None:
            by_child[item.child_id].append(item)
    edge_nodes = {
        node
        for edge in edges
        for node in (edge.left.contig_id, edge.right.contig_id)
    }
    contained: Dict[str, str] = {}
    audits: list[ContainmentAudit] = []
    for child_id in sorted(by_child):
        candidates = by_child[child_id]
        hosts = sorted({item.host_id for item in candidates if item.host_id is not None})
        ranked = sorted(candidates, key=lambda item: (-item.chain.score, item.host_id or ""))
        winner = ranked[0]
        if child_id in edge_nodes:
            reason = "terminal_overlap_conflicts_with_containment"
            status = "rejected"
        elif len(hosts) != 1:
            reason = "competing_containment_hosts"
            status = "rejected"
        else:
            reason = "unique_internal_containment"
            status = "accepted"
            contained[child_id] = hosts[0]
        for item in candidates:
            if item.chain.query_length < item.chain.target_length:
                child_coverage = item.chain.query_coverage
                host_start = item.chain.target_start
                host_end = item.chain.target_end
            else:
                child_coverage = item.chain.target_coverage
                host_start = item.chain.query_start
                host_end = item.chain.query_end
            audits.append(
                ContainmentAudit(
                    child_id=child_id,
                    host_id=item.host_id or ".",
                    strand=item.chain.strand,
                    identity=item.chain.identity,
                    child_coverage=child_coverage,
                    host_start=host_start,
                    host_end=host_end,
                    source_records=tuple(
                        sorted(
                            (record.source, record.line_number)
                            for record in item.chain.records
                        )
                    ),
                    score=item.chain.score,
                    status=status if item is winner else "rejected",
                    reason=reason if item is winner else "duplicate_containment_evidence",
                )
            )
    return contained, tuple(
        sorted(audits, key=lambda row: (row.child_id, row.host_id, -row.score))
    )


def _edge_rank(edge: OrientedOverlapEdge) -> Tuple[object, ...]:
    return (
        -edge.score,
        -min(edge.query_coverage, edge.target_coverage),
        -edge.identity,
        edge.query_id,
        edge.target_id,
        edge.source_records,
    )


def _orient_interval(
    interval: Tuple[int, int],
    length: int,
    orientation: str,
) -> Tuple[int, int]:
    if orientation == "+":
        return interval
    return length - interval[1], length - interval[0]


def _record_signature(record: PafRecord) -> Tuple[object, ...]:
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


def _short_long_coverages(chain: PairChainEvidence) -> Tuple[float, float]:
    if chain.query_length <= chain.target_length:
        return chain.query_coverage, chain.target_coverage
    return chain.target_coverage, chain.query_coverage


def _validate_paf_fasta_boundary(
    records: Sequence[PafRecord],
    sequence_lengths: Mapping[str, int],
) -> None:
    unknown = sorted(
        {
            identifier
            for record in records
            for identifier in (record.query_name, record.target_name)
            if identifier not in sequence_lengths
        }
    )
    if unknown:
        raise ValueError(
            "PAF references IDs absent from FASTA: " + ", ".join(unknown[:5])
        )
    mismatches = sorted(
        {
            record.query_name
            for record in records
            if sequence_lengths[record.query_name] != record.query_length
        }
        | {
            record.target_name
            for record in records
            if sequence_lengths[record.target_name] != record.target_length
        }
    )
    if mismatches:
        raise ValueError(
            "PAF lengths disagree with FASTA for: " + ", ".join(mismatches[:5])
        )
