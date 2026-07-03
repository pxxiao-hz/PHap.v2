"""Build deterministic, oriented alignment-chain evidence from PAF records."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Dict, Tuple

from .locus_rescue import AlignmentChainSummary, interval_union_length
from .paf import PafRecord


@dataclass(frozen=True)
class PafChainEvidence:
    """Auditable evidence for one query-to-target alignment chain."""

    unitig_id: str
    locus_id: str
    query_length: int
    target_length: int
    strand: str
    query_intervals: Tuple[Tuple[int, int], ...]
    target_intervals: Tuple[Tuple[int, int], ...]
    target_start: int
    target_end: int
    matching_bases: int
    alignment_block_bases: int
    identity: float
    union_query_bases: int
    query_coverage: float
    score: float
    collinear: bool
    record_count: int
    records: Tuple[PafRecord, ...]

    def as_locus_summary(self) -> AlignmentChainSummary:
        """Convert to the locus-rescue boundary type."""

        return AlignmentChainSummary(
            unitig_id=self.unitig_id,
            locus_id=self.locus_id,
            query_length=self.query_length,
            query_intervals=self.query_intervals,
            identity=self.identity,
            collinear=self.collinear,
        )


def summarize_paf_chains(
    records: Iterable[PafRecord],
    *,
    max_query_gap: int,
    max_target_gap: int,
) -> Tuple[PafChainEvidence, ...]:
    """Summarize one conservative chain per query/target pair.

    Records for a query/target pair form a supported chain only when they use
    one strand and remain monotonic in both query and oriented target
    coordinates. Gaps are measured in bases. Overlapping intervals are
    retained, but query coverage is calculated from their interval union.
    """

    if max_query_gap < 0 or max_target_gap < 0:
        raise ValueError("chain gap thresholds must be non-negative")

    grouped: Dict[Tuple[str, str], list[PafRecord]] = defaultdict(list)
    query_lengths: Dict[str, int] = {}
    target_lengths: Dict[str, int] = {}
    for record in records:
        previous_query_length = query_lengths.setdefault(
            record.query_name,
            record.query_length,
        )
        if previous_query_length != record.query_length:
            raise ValueError(
                f"inconsistent query length for {record.query_name!r}: "
                f"{previous_query_length} versus {record.query_length}"
            )
        previous_target_length = target_lengths.setdefault(
            record.target_name,
            record.target_length,
        )
        if previous_target_length != record.target_length:
            raise ValueError(
                f"inconsistent target length for {record.target_name!r}: "
                f"{previous_target_length} versus {record.target_length}"
            )
        grouped[(record.query_name, record.target_name)].append(record)

    summaries = [
        _summarize_group(
            unitig_id,
            locus_id,
            group_records,
            max_query_gap=max_query_gap,
            max_target_gap=max_target_gap,
        )
        for (unitig_id, locus_id), group_records in sorted(grouped.items())
    ]
    return tuple(summaries)


def select_longest_targets(
    records: Iterable[PafRecord],
    *,
    target_count: int,
) -> Tuple[str, ...]:
    """Select target IDs by declared length with an ID tie-breaker."""

    if target_count < 1:
        raise ValueError("target_count must be positive")
    target_lengths: Dict[str, int] = {}
    for record in records:
        previous = target_lengths.setdefault(record.target_name, record.target_length)
        if previous != record.target_length:
            raise ValueError(
                f"inconsistent target length for {record.target_name!r}: "
                f"{previous} versus {record.target_length}"
            )
    return tuple(
        target
        for target, _ in sorted(
            target_lengths.items(),
            key=lambda item: (-item[1], item[0]),
        )[:target_count]
    )


def _summarize_group(
    unitig_id: str,
    locus_id: str,
    records: Iterable[PafRecord],
    *,
    max_query_gap: int,
    max_target_gap: int,
) -> PafChainEvidence:
    ordered = tuple(sorted(records, key=_record_sort_key))
    if not ordered:
        raise ValueError("cannot summarize an empty PAF chain")
    query_length = ordered[0].query_length
    target_length = ordered[0].target_length
    strands = {record.strand for record in ordered}
    strand = next(iter(strands)) if len(strands) == 1 else "mixed"
    collinear = len(strands) == 1 and _is_collinear(
        ordered,
        strand,
        max_query_gap=max_query_gap,
        max_target_gap=max_target_gap,
    )
    query_intervals = tuple(
        (record.query_start, record.query_end) for record in ordered
    )
    target_intervals = tuple(
        (record.target_start, record.target_end) for record in ordered
    )
    union_query_bases = interval_union_length(
        query_intervals,
        sequence_length=query_length,
    )
    block_bases = sum(record.alignment_block_length for record in ordered)
    matching_bases = sum(record.matching_bases for record in ordered)
    identity = matching_bases / block_bases if block_bases else 0.0
    if not math.isfinite(identity) or not 0.0 <= identity <= 1.0:
        raise AssertionError("validated PAF records produced an invalid identity")
    query_coverage = union_query_bases / query_length
    return PafChainEvidence(
        unitig_id=unitig_id,
        locus_id=locus_id,
        query_length=query_length,
        target_length=target_length,
        strand=strand,
        query_intervals=query_intervals,
        target_intervals=target_intervals,
        target_start=min(record.target_start for record in ordered),
        target_end=max(record.target_end for record in ordered),
        matching_bases=matching_bases,
        alignment_block_bases=block_bases,
        identity=identity,
        union_query_bases=union_query_bases,
        query_coverage=query_coverage,
        score=identity * query_coverage,
        collinear=collinear,
        record_count=len(ordered),
        records=ordered,
    )


def _is_collinear(
    records: Tuple[PafRecord, ...],
    strand: str,
    *,
    max_query_gap: int,
    max_target_gap: int,
) -> bool:
    for previous, current in zip(records, records[1:]):
        query_gap = current.query_start - previous.query_end
        if query_gap > max_query_gap:
            return False
        if strand == "+":
            if current.target_start < previous.target_start:
                return False
            target_gap = current.target_start - previous.target_end
        elif strand == "-":
            if current.target_start > previous.target_start:
                return False
            target_gap = previous.target_start - current.target_end
        else:
            return False
        if target_gap > max_target_gap:
            return False
    return True


def _record_sort_key(record: PafRecord) -> Tuple[object, ...]:
    return (
        record.query_start,
        record.query_end,
        record.target_start,
        record.target_end,
        record.strand,
        record.line_number,
    )
