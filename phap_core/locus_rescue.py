"""Pure, auditable placement of unsupported unitigs on an mT2T reference.

This module consumes alignment-chain summaries produced elsewhere.  It does
not parse PAF or inspect BAM files.  Query intervals must use 0-based,
half-open coordinates.

An mT2T reference can establish a chromosome/bin (``locus_id``), but because
the reference is a mixed-haplotype sequence it cannot establish a haplotype
group by itself.  A group is returned only when independent, group-specific
read or overlap evidence supports exactly one compatible group.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Literal, Optional, Sequence, Tuple


RescueStatus = Literal[
    "locus_assigned_haplotype_unresolved",
    "rescued_group",
    "ambiguous",
    "unassigned",
]


class LocusRescueError(ValueError):
    """Raised when locus-rescue evidence or thresholds are invalid."""


@dataclass(frozen=True)
class AlignmentChainSummary:
    """One precomputed unitig-to-locus alignment chain.

    ``identity`` is a fraction in ``[0, 1]``. ``query_intervals`` may overlap;
    query coverage is calculated from their interval union rather than their
    summed lengths.
    """

    unitig_id: str
    locus_id: str
    query_length: int
    query_intervals: Tuple[Tuple[int, int], ...]
    identity: float
    collinear: bool = True


@dataclass(frozen=True)
class GroupSpecificEvidence:
    """Independent evidence connecting a unitig to one haplotype group."""

    group_id: str
    locus_id: str
    read_support: int = 0
    overlap_support: float = 0.0


@dataclass(frozen=True)
class LocusRescueThresholds:
    """Explicit thresholds used by :func:`decide_locus_rescue`.

    Identity, query coverage, and the next-best score margin are fractions.
    The score used for locus ranking is ``identity * query_coverage``.
    """

    min_identity: float
    min_query_coverage: float
    min_next_best_margin: float
    min_low_coverage_read_support: int = 1
    min_group_read_support: int = 1
    min_group_overlap_support: float = 1.0

    def __post_init__(self) -> None:
        for name, value in (
            ("min_identity", self.min_identity),
            ("min_query_coverage", self.min_query_coverage),
            ("min_next_best_margin", self.min_next_best_margin),
        ):
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise LocusRescueError(f"{name} must be a finite fraction in [0, 1]")
        if self.min_low_coverage_read_support < 0:
            raise LocusRescueError("min_low_coverage_read_support must be non-negative")
        if self.min_group_read_support < 0:
            raise LocusRescueError("min_group_read_support must be non-negative")
        if (
            not math.isfinite(self.min_group_overlap_support)
            or self.min_group_overlap_support < 0
        ):
            raise LocusRescueError(
                "min_group_overlap_support must be finite and non-negative"
            )


@dataclass(frozen=True)
class LocusCandidateAudit:
    """Calculated evidence for one candidate mT2T locus."""

    locus_id: str
    identity: float
    query_coverage: float
    union_query_bases: int
    score: float
    collinear: bool
    meets_thresholds: bool


@dataclass(frozen=True)
class LocusRescueDecision:
    """One typed, tabular-ready rescue decision."""

    unitig_id: str
    source_state: str
    status: RescueStatus
    assigned_locus: Optional[str]
    assigned_group: Optional[str]
    rescue_class: Optional[str]
    reason: str
    independent_read_support: int
    best_locus: Optional[str]
    best_identity: Optional[float]
    best_query_coverage: Optional[float]
    best_score: Optional[float]
    next_best_locus: Optional[str]
    next_best_score: Optional[float]
    next_best_margin: Optional[float]
    supported_groups: Tuple[str, ...]
    candidates: Tuple[LocusCandidateAudit, ...]


def interval_union_length(
    intervals: Iterable[Tuple[int, int]],
    *,
    sequence_length: Optional[int] = None,
) -> int:
    """Return covered bases after merging 0-based half-open intervals."""

    if sequence_length is not None and sequence_length < 1:
        raise LocusRescueError("sequence_length must be positive")

    validated: list[Tuple[int, int]] = []
    for start, end in intervals:
        if start < 0 or end < start:
            raise LocusRescueError(f"invalid half-open interval: ({start}, {end})")
        if sequence_length is not None and end > sequence_length:
            raise LocusRescueError(
                f"interval ({start}, {end}) exceeds sequence length {sequence_length}"
            )
        if start != end:
            validated.append((start, end))
    if not validated:
        return 0

    total = 0
    current_start, current_end = sorted(validated)[0]
    for start, end in sorted(validated)[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    return total + current_end - current_start


def summarize_locus_candidates(
    unitig_id: str,
    alignments: Sequence[AlignmentChainSummary],
    thresholds: LocusRescueThresholds,
) -> Tuple[LocusCandidateAudit, ...]:
    """Validate and rank precomputed alignment chains deterministically."""

    if not unitig_id:
        raise LocusRescueError("unitig_id must not be empty")

    seen_loci: set[str] = set()
    audits: list[LocusCandidateAudit] = []
    for alignment in alignments:
        if alignment.unitig_id != unitig_id:
            raise LocusRescueError(
                f"alignment unitig {alignment.unitig_id!r} does not match {unitig_id!r}"
            )
        if not alignment.locus_id:
            raise LocusRescueError("locus_id must not be empty")
        if alignment.locus_id in seen_loci:
            raise LocusRescueError(
                f"duplicate alignment summary for locus {alignment.locus_id!r}"
            )
        seen_loci.add(alignment.locus_id)
        if alignment.query_length < 1:
            raise LocusRescueError("query_length must be positive")
        if not math.isfinite(alignment.identity) or not 0.0 <= alignment.identity <= 1.0:
            raise LocusRescueError("alignment identity must be a finite fraction in [0, 1]")

        union_bases = interval_union_length(
            alignment.query_intervals,
            sequence_length=alignment.query_length,
        )
        query_coverage = union_bases / alignment.query_length
        score = alignment.identity * query_coverage
        meets_thresholds = (
            alignment.collinear
            and alignment.identity >= thresholds.min_identity
            and query_coverage >= thresholds.min_query_coverage
        )
        audits.append(
            LocusCandidateAudit(
                locus_id=alignment.locus_id,
                identity=alignment.identity,
                query_coverage=query_coverage,
                union_query_bases=union_bases,
                score=score,
                collinear=alignment.collinear,
                meets_thresholds=meets_thresholds,
            )
        )

    return tuple(
        sorted(
            audits,
            key=lambda audit: (
                -audit.score,
                -audit.query_coverage,
                -audit.identity,
                audit.locus_id,
            ),
        )
    )


def decide_locus_rescue(
    unitig_id: str,
    alignments: Sequence[AlignmentChainSummary],
    *,
    source_state: str,
    independent_read_support: int,
    thresholds: LocusRescueThresholds,
    group_evidence: Sequence[GroupSpecificEvidence] = (),
) -> LocusRescueDecision:
    """Place a low-support unitig without inferring haplotype from mT2T alone."""

    if not source_state:
        raise LocusRescueError("source_state must not be empty")
    if independent_read_support < 0:
        raise LocusRescueError("independent_read_support must be non-negative")

    candidates = summarize_locus_candidates(unitig_id, alignments, thresholds)
    eligible = tuple(candidate for candidate in candidates if candidate.meets_thresholds)
    if not eligible:
        return _decision(
            unitig_id,
            source_state,
            "unassigned",
            "no_collinear_locus_meets_thresholds",
            independent_read_support,
            candidates,
        )

    best = eligible[0]
    competitors = tuple(
        candidate
        for candidate in candidates
        if candidate.collinear and candidate.locus_id != best.locus_id
    )
    next_best = competitors[0] if competitors else None
    margin = best.score - next_best.score if next_best is not None else None
    if margin is not None and margin < thresholds.min_next_best_margin:
        return _decision(
            unitig_id,
            source_state,
            "ambiguous",
            "competing_mT2T_loci",
            independent_read_support,
            candidates,
            best=best,
            next_best=next_best,
            margin=margin,
        )

    low_coverage = source_state in {
        "low_coverage",
        "low_coverage_haplotig_candidate",
    }
    if (
        low_coverage
        and independent_read_support < thresholds.min_low_coverage_read_support
    ):
        return _decision(
            unitig_id,
            source_state,
            "unassigned",
            "low_coverage_without_independent_read_support",
            independent_read_support,
            candidates,
            best=best,
            next_best=next_best,
            margin=margin,
        )

    supported = _supported_group_evidence(group_evidence, thresholds)
    conflicting_loci = tuple(
        evidence for evidence in supported if evidence.locus_id != best.locus_id
    )
    compatible_groups = tuple(
        sorted(
            {
                evidence.group_id
                for evidence in supported
                if evidence.locus_id == best.locus_id
            }
        )
    )
    rescue_class = "low_coverage_haplotig_candidate" if low_coverage else None
    if conflicting_loci or len(compatible_groups) > 1:
        return _decision(
            unitig_id,
            source_state,
            "ambiguous",
            "conflicting_group_specific_evidence",
            independent_read_support,
            candidates,
            best=best,
            next_best=next_best,
            margin=margin,
            rescue_class=rescue_class,
            supported_groups=compatible_groups,
        )
    if len(compatible_groups) == 1:
        return _decision(
            unitig_id,
            source_state,
            "rescued_group",
            "unique_locus_with_group_specific_evidence",
            independent_read_support,
            candidates,
            best=best,
            next_best=next_best,
            margin=margin,
            assigned_locus=best.locus_id,
            assigned_group=compatible_groups[0],
            rescue_class=rescue_class,
            supported_groups=compatible_groups,
        )

    return _decision(
        unitig_id,
        source_state,
        "locus_assigned_haplotype_unresolved",
        "mT2T_supports_locus_not_haplotype",
        independent_read_support,
        candidates,
        best=best,
        next_best=next_best,
        margin=margin,
        assigned_locus=best.locus_id,
        rescue_class=rescue_class,
    )


def _supported_group_evidence(
    evidence_rows: Sequence[GroupSpecificEvidence],
    thresholds: LocusRescueThresholds,
) -> Tuple[GroupSpecificEvidence, ...]:
    supported: list[GroupSpecificEvidence] = []
    seen_pairs: set[Tuple[str, str]] = set()
    for evidence in evidence_rows:
        if not evidence.group_id or not evidence.locus_id:
            raise LocusRescueError("group and locus identifiers must not be empty")
        pair = (evidence.group_id, evidence.locus_id)
        if pair in seen_pairs:
            raise LocusRescueError(f"duplicate group-specific evidence for {pair!r}")
        seen_pairs.add(pair)
        if evidence.read_support < 0:
            raise LocusRescueError("group read support must be non-negative")
        if not math.isfinite(evidence.overlap_support) or evidence.overlap_support < 0:
            raise LocusRescueError("group overlap support must be finite and non-negative")
        if (
            evidence.read_support >= thresholds.min_group_read_support
            or evidence.overlap_support >= thresholds.min_group_overlap_support
        ):
            supported.append(evidence)
    return tuple(sorted(supported, key=lambda row: (row.locus_id, row.group_id)))


def _decision(
    unitig_id: str,
    source_state: str,
    status: RescueStatus,
    reason: str,
    independent_read_support: int,
    candidates: Tuple[LocusCandidateAudit, ...],
    *,
    best: Optional[LocusCandidateAudit] = None,
    next_best: Optional[LocusCandidateAudit] = None,
    margin: Optional[float] = None,
    assigned_locus: Optional[str] = None,
    assigned_group: Optional[str] = None,
    rescue_class: Optional[str] = None,
    supported_groups: Tuple[str, ...] = (),
) -> LocusRescueDecision:
    return LocusRescueDecision(
        unitig_id=unitig_id,
        source_state=source_state,
        status=status,
        assigned_locus=assigned_locus,
        assigned_group=assigned_group,
        rescue_class=rescue_class,
        reason=reason,
        independent_read_support=independent_read_support,
        best_locus=best.locus_id if best is not None else None,
        best_identity=best.identity if best is not None else None,
        best_query_coverage=best.query_coverage if best is not None else None,
        best_score=best.score if best is not None else None,
        next_best_locus=next_best.locus_id if next_best is not None else None,
        next_best_score=next_best.score if next_best is not None else None,
        next_best_margin=margin,
        supported_groups=supported_groups,
        candidates=candidates,
    )
