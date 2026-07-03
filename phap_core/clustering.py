"""Deterministic dosage-aware clustering over an mT2T allelic table."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Optional, Tuple

from .hic_scoring import (
    HiCGroupScore,
    ScoreMode,
    assign_candidate_groups,
    validate_pair_links,
)


@dataclass(frozen=True)
class AllelicBin:
    """One mT2T interval and the unitigs aligned to it."""

    chromosome: str
    start: int
    end: int
    unitigs: Tuple[str, ...]


@dataclass(frozen=True)
class ClusterDecision:
    """Final auditable destination of one unitig within one chromosome."""

    unitig_id: str
    chromosome: str
    start: int
    end: int
    source_state: str
    dosage: Optional[int]
    status: str
    selected_groups: Tuple[str, ...]
    reason: str
    score_mode: ScoreMode
    scores: Tuple[HiCGroupScore, ...]


@dataclass(frozen=True)
class ChromosomeClusterResult:
    """Groups and one final decision for every observed unitig."""

    group_unitigs: Tuple[Tuple[str, Tuple[str, ...]], ...]
    decisions: Tuple[ClusterDecision, ...]


def cluster_allelic_bins(
    rows: Sequence[AllelicBin],
    dosage_by_unitig: Mapping[str, Optional[int]],
    source_states: Mapping[str, str],
    pair_links: Mapping[Tuple[str, str], float],
    restriction_sites: Mapping[str, int],
    *,
    ploidy: int,
    score_mode: ScoreMode,
    min_score: float,
    min_margin: float,
) -> ChromosomeClusterResult:
    """Cluster unitigs while conserving dosage and retaining unsupported records.

    The first valid allelic bin seeds anonymous haplotype groups
    deterministically.  Later unitigs use Hi-C evidence against the available
    groups in their mT2T bin.  Existing unitig membership supplies path
    continuity.  A failed decision never produces a partial assignment.
    """

    if ploidy < 1:
        raise ValueError("ploidy must be positive")
    if not rows:
        raise ValueError("at least one allelic bin is required")
    chromosomes = {row.chromosome for row in rows}
    if len(chromosomes) != 1 or "" in chromosomes:
        raise ValueError("rows must contain exactly one non-empty chromosome")
    _validate_rows(rows)
    validate_pair_links(pair_links)

    group_ids = tuple(f"group{index}" for index in range(1, ploidy + 1))
    groups: dict[str, set[str]] = {group_id: set() for group_id in group_ids}
    memberships: dict[str, Tuple[str, ...]] = {}
    decisions: dict[str, ClusterDecision] = {}

    for row in rows:
        seeding_row = not any(groups.values())
        occupied = {
            group_id
            for unitig_id in row.unitigs
            for group_id in memberships.get(unitig_id, ())
        }
        new_unitigs = sorted(
            {unitig_id for unitig_id in row.unitigs if unitig_id not in memberships},
            key=lambda unitig_id: (
                -_sortable_dosage(dosage_by_unitig.get(unitig_id)),
                unitig_id,
            ),
        )
        for unitig_id in new_unitigs:
            dosage = dosage_by_unitig.get(unitig_id)
            source_state = source_states.get(unitig_id, "missing_dosage")
            if dosage is None or not 1 <= dosage <= ploidy:
                decisions[unitig_id] = _unresolved_decision(
                    unitig_id,
                    row,
                    source_state,
                    dosage,
                    score_mode,
                    "invalid_or_missing_dosage",
                )
                continue

            available = tuple(
                group_id for group_id in group_ids if group_id not in occupied
            )
            if len(available) < dosage:
                decisions[unitig_id] = _unresolved_decision(
                    unitig_id,
                    row,
                    source_state,
                    dosage,
                    score_mode,
                    "insufficient_bin_capacity",
                )
                continue

            if seeding_row:
                selected = available[:dosage]
                decision = ClusterDecision(
                    unitig_id=unitig_id,
                    chromosome=row.chromosome,
                    start=row.start,
                    end=row.end,
                    source_state=source_state,
                    dosage=dosage,
                    status="assigned",
                    selected_groups=selected,
                    reason="deterministic_initial_seed",
                    score_mode=score_mode,
                    scores=(),
                )
            else:
                candidate_groups = {
                    group_id: groups[group_id] for group_id in available
                }
                links_by_unitig = _links_for_candidate(unitig_id, pair_links)
                assignment = assign_candidate_groups(
                    unitig_id,
                    candidate_groups,
                    links_by_unitig,
                    restriction_sites,
                    dosage=dosage,
                    ploidy=len(candidate_groups),
                    mode=score_mode,
                    min_score=min_score,
                    min_margin=min_margin,
                )
                if assignment.status != "assigned":
                    raw_support = sum(score.raw_score for score in assignment.scores)
                    reason = (
                        "no_hic_support"
                        if raw_support == 0
                        else assignment.reason
                    )
                    decisions[unitig_id] = _unresolved_decision(
                        unitig_id,
                        row,
                        source_state,
                        dosage,
                        score_mode,
                        reason,
                        scores=assignment.scores,
                    )
                    continue
                selected = assignment.selected_groups
                decision = ClusterDecision(
                    unitig_id=unitig_id,
                    chromosome=row.chromosome,
                    start=row.start,
                    end=row.end,
                    source_state=source_state,
                    dosage=dosage,
                    status="assigned",
                    selected_groups=selected,
                    reason=assignment.reason,
                    score_mode=score_mode,
                    scores=assignment.scores,
                )

            if len(selected) != dosage:
                raise AssertionError("assigned unitig does not conserve dosage")
            for group_id in selected:
                groups[group_id].add(unitig_id)
            memberships[unitig_id] = tuple(sorted(selected))
            occupied.update(selected)
            decisions[unitig_id] = decision

    observed = {unitig_id for row in rows for unitig_id in row.unitigs}
    if observed != set(decisions):
        raise AssertionError("every observed unitig must have a final decision")
    return ChromosomeClusterResult(
        group_unitigs=tuple(
            (group_id, tuple(sorted(groups[group_id]))) for group_id in group_ids
        ),
        decisions=tuple(decisions[unitig_id] for unitig_id in sorted(decisions)),
    )


def _unresolved_decision(
    unitig_id: str,
    row: AllelicBin,
    source_state: str,
    dosage: Optional[int],
    score_mode: ScoreMode,
    reason: str,
    *,
    scores: Tuple[HiCGroupScore, ...] = (),
) -> ClusterDecision:
    status = (
        "locus_assigned_haplotype_unresolved"
        if reason in {"no_hic_support", "invalid_or_missing_dosage"}
        else "ambiguous"
    )
    return ClusterDecision(
        unitig_id=unitig_id,
        chromosome=row.chromosome,
        start=row.start,
        end=row.end,
        source_state=source_state,
        dosage=dosage,
        status=status,
        selected_groups=(),
        reason=reason,
        score_mode=score_mode,
        scores=scores,
    )


def _links_for_candidate(
    unitig_id: str,
    pair_links: Mapping[Tuple[str, str], float],
) -> dict[str, float]:
    links: dict[str, float] = {}
    for (left, right), count in pair_links.items():
        if left == unitig_id and right != unitig_id:
            links[right] = links.get(right, 0.0) + count
        elif right == unitig_id and left != unitig_id:
            links[left] = links.get(left, 0.0) + count
    return links


def _sortable_dosage(dosage: Optional[int]) -> int:
    return dosage if dosage is not None and dosage > 0 else 0


def _validate_rows(rows: Sequence[AllelicBin]) -> None:
    for row in rows:
        if row.start < 0 or row.end < row.start:
            raise ValueError("allelic-bin coordinates must be 0-based half-open")
        if not row.unitigs or any(not unitig_id for unitig_id in row.unitigs):
            raise ValueError("each allelic bin must contain non-empty unitig IDs")
