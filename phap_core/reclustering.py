"""Pure reassignment of unclustered unitigs with chromosome-aware Hi-C scores."""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Optional, Tuple

from .hic_scoring import (
    HiCGroupAssignment,
    HiCGroupScore,
    ScoreMode,
    assign_candidate_groups,
    validate_pair_links,
)


_GROUP_PATTERN = re.compile(r"^(?:(?P<family>.+)_)?group(?P<index>[1-9][0-9]*)$")


@dataclass(frozen=True)
class ReclusterDecision:
    unitig_id: str
    source_state: str
    dosage: Optional[int]
    status: str
    selected_groups: Tuple[str, ...]
    selected_family: Optional[str]
    assigned_locus: Optional[str]
    reason: str
    scores: Tuple[HiCGroupScore, ...]


@dataclass(frozen=True)
class ReclusterResult:
    memberships: Tuple[Tuple[str, Tuple[str, ...]], ...]
    decisions: Tuple[ReclusterDecision, ...]


def reassign_unitigs(
    unitig_ids: Sequence[str],
    existing_memberships: Mapping[str, Iterable[str]],
    dosage_by_unitig: Mapping[str, Optional[int]],
    source_states: Mapping[str, str],
    pair_links: Mapping[Tuple[str, str], float],
    restriction_sites: Mapping[str, int],
    *,
    ploidy: int,
    score_mode: ScoreMode,
    min_score: float,
    min_margin: float,
    known_locus: Optional[str] = None,
    declared_groups: Optional[Sequence[str]] = None,
) -> ReclusterResult:
    """Assign exactly ``dosage`` groups within one chromosome family or none."""

    if ploidy < 1:
        raise ValueError("ploidy must be positive")
    validate_pair_links(pair_links)
    memberships = {
        unitig_id: set(groups) for unitig_id, groups in existing_memberships.items()
    }
    _validate_existing_memberships(memberships, dosage_by_unitig, ploidy)
    group_members = _invert_memberships(memberships, declared_groups)
    families = _group_families(group_members, ploidy)
    decisions: list[ReclusterDecision] = []

    for unitig_id in unitig_ids:
        if not unitig_id:
            raise ValueError("unitig IDs must not be empty")
        source_state = source_states.get(unitig_id, "missing_dosage")
        dosage = dosage_by_unitig.get(unitig_id)
        if unitig_id in memberships:
            decisions.append(
                ReclusterDecision(
                    unitig_id,
                    source_state,
                    dosage,
                    "assigned",
                    tuple(sorted(memberships[unitig_id])),
                    _single_membership_family(memberships[unitig_id]),
                    known_locus,
                    "preassigned",
                    (),
                )
            )
            continue
        if dosage is None or not 1 <= dosage <= ploidy:
            decisions.append(
                _unresolved(
                    unitig_id,
                    source_state,
                    dosage,
                    "invalid_or_missing_dosage",
                    known_locus,
                )
            )
            continue

        links = _links_for_candidate(unitig_id, pair_links)
        family_assignments = [
            (
                family,
                assign_candidate_groups(
                    unitig_id,
                    {group: group_members[group] for group in groups},
                    links,
                    restriction_sites,
                    dosage=dosage,
                    ploidy=ploidy,
                    mode=score_mode,
                    min_score=min_score,
                    min_margin=min_margin,
                ),
            )
            for family, groups in families
        ]
        all_scores = tuple(
            score
            for _, assignment in family_assignments
            for score in assignment.scores
        )
        supported = [
            (family, assignment, _selected_score(assignment))
            for family, assignment in family_assignments
            if assignment.status == "assigned"
        ]
        if not supported:
            raw_support = sum(score.raw_score for score in all_scores)
            reason = "no_hic_support" if raw_support == 0 else "no_supported_group_family"
            decisions.append(
                _unresolved(
                    unitig_id,
                    source_state,
                    dosage,
                    reason,
                    known_locus,
                    scores=all_scores,
                )
            )
            continue

        supported.sort(key=lambda item: (-item[2], item[0]))
        best_family, best_assignment, best_score = supported[0]
        if len(supported) > 1:
            family_margin = best_score - supported[1][2]
            if family_margin == 0 or family_margin < min_margin:
                decisions.append(
                    ReclusterDecision(
                        unitig_id,
                        source_state,
                        dosage,
                        "ambiguous",
                        (),
                        None,
                        known_locus,
                        "competing_chromosome_families",
                        all_scores,
                    )
                )
                continue

        selected = best_assignment.selected_groups
        if len(selected) != dosage:
            raise AssertionError("recluster assignment does not conserve dosage")
        memberships[unitig_id] = set(selected)
        for group_id in selected:
            group_members[group_id].add(unitig_id)
        decisions.append(
            ReclusterDecision(
                unitig_id,
                source_state,
                dosage,
                "assigned",
                selected,
                best_family or None,
                known_locus or best_family or None,
                best_assignment.reason,
                all_scores,
            )
        )

    return ReclusterResult(
        memberships=tuple(
            (unitig_id, tuple(sorted(groups)))
            for unitig_id, groups in sorted(memberships.items())
        ),
        decisions=tuple(decisions),
    )


def _invert_memberships(
    memberships: Mapping[str, set[str]],
    declared_groups: Optional[Sequence[str]],
) -> dict[str, set[str]]:
    groups: dict[str, set[str]] = {
        group_id: set() for group_id in (declared_groups or ())
    }
    if any(not group_id for group_id in groups):
        raise ValueError("declared group IDs must not be empty")
    for unitig_id, group_ids in memberships.items():
        for group_id in group_ids:
            groups.setdefault(group_id, set()).add(unitig_id)
    if not groups:
        raise ValueError("at least one pre-existing group is required")
    return groups


def _validate_existing_memberships(
    memberships: Mapping[str, set[str]],
    dosage_by_unitig: Mapping[str, Optional[int]],
    ploidy: int,
) -> None:
    for unitig_id, groups in memberships.items():
        dosage = dosage_by_unitig.get(unitig_id)
        if dosage is None or not 1 <= dosage <= ploidy:
            raise ValueError(
                f"preassigned unitig {unitig_id!r} has invalid or missing dosage"
            )
        if len(groups) != dosage:
            raise ValueError(
                f"preassigned unitig {unitig_id!r} occupies {len(groups)} groups, "
                f"expected dosage {dosage}"
            )


def _group_families(
    group_members: Mapping[str, set[str]],
    ploidy: int,
) -> Tuple[Tuple[str, Tuple[str, ...]], ...]:
    by_family: dict[str, list[Tuple[int, str]]] = {}
    for group_id in group_members:
        match = _GROUP_PATTERN.fullmatch(group_id)
        if match is None:
            raise ValueError(f"unsupported group ID format: {group_id!r}")
        family = match.group("family") or ""
        by_family.setdefault(family, []).append((int(match.group("index")), group_id))
    result = []
    for family in sorted(by_family):
        indexed = sorted(by_family[family])
        if [index for index, _ in indexed] != list(range(1, ploidy + 1)):
            raise ValueError(
                f"group family {family or '<default>'!r} must contain group1..group{ploidy}"
            )
        result.append((family, tuple(group_id for _, group_id in indexed)))
    return tuple(result)


def _single_membership_family(groups: set[str]) -> Optional[str]:
    families = {
        (match.group("family") or "")
        for group_id in groups
        if (match := _GROUP_PATTERN.fullmatch(group_id)) is not None
    }
    return next(iter(families)) if len(families) == 1 else None


def _selected_score(assignment: HiCGroupAssignment) -> float:
    selected = set(assignment.selected_groups)
    if assignment.mode == "raw":
        return sum(
            score.raw_score for score in assignment.scores if score.group_id in selected
        )
    return sum(
        score.re_density or 0.0
        for score in assignment.scores
        if score.group_id in selected
    )


def _unresolved(
    unitig_id: str,
    source_state: str,
    dosage: Optional[int],
    reason: str,
    known_locus: Optional[str],
    *,
    scores: Tuple[HiCGroupScore, ...] = (),
) -> ReclusterDecision:
    if known_locus is not None and reason in {
        "invalid_or_missing_dosage",
        "no_hic_support",
    }:
        status = "locus_assigned_haplotype_unresolved"
    elif reason in {"invalid_or_missing_dosage", "no_hic_support"}:
        status = "unassigned"
    else:
        status = "ambiguous"
    return ReclusterDecision(
        unitig_id,
        source_state,
        dosage,
        status,
        (),
        None,
        known_locus,
        reason,
        scores,
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
