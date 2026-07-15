"""Dosage-capacity selection for mT2T allelic-bin candidates."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import DefaultDict, Dict, FrozenSet, Literal, Optional, Tuple

from .allelic_bins import AllelicBinCandidate, AllelicBinKey


CandidateStatus = Literal["selected", "rejected", "ambiguous", "unassigned"]
BinStatus = Literal["selected", "partial_capacity", "ambiguous", "unassigned"]


@dataclass(frozen=True)
class AllelicCandidateDecision:
    """One auditable routing decision for one bin/unitig candidate."""

    candidate: AllelicBinCandidate
    source_state: str
    dosage: Optional[int]
    copy_weighted_support_bases: Optional[int]
    status: CandidateStatus
    reason: str


@dataclass(frozen=True)
class AllelicBinDecision:
    """Capacity-constrained decision for one target bin."""

    bin_key: AllelicBinKey
    ploidy: int
    min_support_bases: int
    min_bin_coverage: float
    status: BinStatus
    reason: str
    selected_unitigs: Tuple[str, ...]
    selected_copy_count: int
    remaining_capacity: int
    best_objective: Optional[int]
    optimal_solution_count: int
    candidates: Tuple[AllelicCandidateDecision, ...]


@dataclass(frozen=True)
class _SolutionState:
    score: int
    count: int
    representative: Tuple[str, ...]
    member_union: FrozenSet[str]


def select_dosage_capacity_bins(
    candidates: Sequence[AllelicBinCandidate],
    dosage_by_unitig: Mapping[str, Optional[int]],
    source_states: Mapping[str, str],
    *,
    ploidy: int,
    min_support_bases: int = 1,
    min_bin_coverage: float = 0.0,
) -> Tuple[AllelicBinDecision, ...]:
    """Select a unique best subset per bin under the copy-count capacity.

    The integer objective is the sum of ``target_union_bases * dosage``.  It
    compares a collapsed unitig with the same number of independently assembled
    haplotigs in copy-base units.  Equal best subsets are reported as ambiguous;
    identifiers are used only to stabilize output order, never to choose one.
    """

    if ploidy < 1:
        raise ValueError("ploidy must be positive")
    if min_support_bases < 1:
        raise ValueError("min_support_bases must be positive")
    if not math.isfinite(min_bin_coverage) or not 0.0 <= min_bin_coverage <= 1.0:
        raise ValueError("min_bin_coverage must be a finite fraction in [0, 1]")

    by_bin: DefaultDict[AllelicBinKey, list[AllelicBinCandidate]] = defaultdict(list)
    seen: set[Tuple[AllelicBinKey, str]] = set()
    for candidate in candidates:
        key = candidate.bin_key
        if not key.target_id or key.start < 0 or key.end <= key.start:
            raise ValueError("invalid allelic bin key")
        if not candidate.unitig_id:
            raise ValueError("unitig_id must not be empty")
        if not 0 < candidate.union_support_bases <= key.end - key.start:
            raise ValueError("union support must be within the bin length")
        identity = (candidate.bin_key, candidate.unitig_id)
        if identity in seen:
            raise ValueError(
                "duplicate allelic-bin candidate for "
                f"{candidate.bin_key.target_id}:{candidate.bin_key.start}-"
                f"{candidate.bin_key.end}:{candidate.unitig_id}"
            )
        seen.add(identity)
        by_bin[candidate.bin_key].append(candidate)

    return tuple(
        _select_one_bin(
            key,
            by_bin[key],
            dosage_by_unitig,
            source_states,
            ploidy=ploidy,
            min_support_bases=min_support_bases,
            min_bin_coverage=min_bin_coverage,
        )
        for key in sorted(by_bin)
    )


def _select_one_bin(
    key: AllelicBinKey,
    candidates: Sequence[AllelicBinCandidate],
    dosage_by_unitig: Mapping[str, Optional[int]],
    source_states: Mapping[str, str],
    *,
    ploidy: int,
    min_support_bases: int,
    min_bin_coverage: float,
) -> AllelicBinDecision:
    ordered = tuple(sorted(candidates, key=lambda candidate: candidate.unitig_id))
    eligible: list[Tuple[AllelicBinCandidate, int, int]] = []
    invalid: Dict[str, Tuple[Optional[int], str]] = {}
    for candidate in ordered:
        dosage = dosage_by_unitig.get(candidate.unitig_id)
        if candidate.unitig_id not in dosage_by_unitig or dosage is None:
            invalid[candidate.unitig_id] = (dosage, "invalid_or_missing_dosage")
        elif not 1 <= dosage <= ploidy:
            invalid[candidate.unitig_id] = (dosage, "dosage_outside_ploidy")
        elif candidate.union_support_bases < min_support_bases:
            invalid[candidate.unitig_id] = (dosage, "insufficient_support_bases")
        elif (
            candidate.union_support_bases / (key.end - key.start)
            < min_bin_coverage
        ):
            invalid[candidate.unitig_id] = (dosage, "insufficient_bin_coverage")
        else:
            eligible.append(
                (
                    candidate,
                    dosage,
                    candidate.union_support_bases * dosage,
                )
            )

    if not eligible:
        invalid_decisions = tuple(
            AllelicCandidateDecision(
                candidate=candidate,
                source_state=source_states.get(candidate.unitig_id, "missing_dosage"),
                dosage=invalid[candidate.unitig_id][0],
                copy_weighted_support_bases=None,
                status="unassigned",
                reason=invalid[candidate.unitig_id][1],
            )
            for candidate in ordered
        )
        return AllelicBinDecision(
            bin_key=key,
            ploidy=ploidy,
            min_support_bases=min_support_bases,
            min_bin_coverage=min_bin_coverage,
            status="unassigned",
            reason="no_eligible_candidate",
            selected_unitigs=(),
            selected_copy_count=0,
            remaining_capacity=ploidy,
            best_objective=None,
            optimal_solution_count=0,
            candidates=invalid_decisions,
        )

    states: Dict[int, _SolutionState] = {
        0: _SolutionState(0, 1, (), frozenset())
    }
    eligible_by_id = {
        candidate.unitig_id: (candidate, dosage, score)
        for candidate, dosage, score in eligible
    }
    for unitig_id in sorted(eligible_by_id):
        _, dosage, score = eligible_by_id[unitig_id]
        updated = dict(states)
        for used_capacity, state in states.items():
            new_capacity = used_capacity + dosage
            if new_capacity > ploidy:
                continue
            proposal = _SolutionState(
                score=state.score + score,
                count=state.count,
                representative=tuple(sorted((*state.representative, unitig_id))),
                member_union=state.member_union | {unitig_id},
            )
            incumbent = updated.get(new_capacity)
            if incumbent is None or proposal.score > incumbent.score:
                updated[new_capacity] = proposal
            elif proposal.score == incumbent.score:
                updated[new_capacity] = _merge_equal_states(incumbent, proposal)
        states = updated

    nonempty_states = [
        (capacity, state) for capacity, state in states.items() if capacity > 0
    ]
    best_score = max(state.score for _, state in nonempty_states)
    optimal = [
        (capacity, state)
        for capacity, state in nonempty_states
        if state.score == best_score
    ]
    optimal_count = sum(state.count for _, state in optimal)
    optimal_members = frozenset().union(
        *(state.member_union for _, state in optimal)
    )

    if optimal_count > 1:
        selected_unitigs: Tuple[str, ...] = ()
        selected_copy_count = 0
        status: BinStatus = "ambiguous"
        reason = "multiple_equal_optimum_subsets"
    else:
        selected_capacity, selected_state = optimal[0]
        selected_unitigs = selected_state.representative
        selected_copy_count = selected_capacity
        status = "selected" if selected_capacity == ploidy else "partial_capacity"
        reason = (
            "unique_optimum_fills_ploidy"
            if selected_capacity == ploidy
            else "unique_optimum_partial_capacity"
        )

    candidate_decisions: list[AllelicCandidateDecision] = []
    for candidate in ordered:
        unitig_id = candidate.unitig_id
        source_state = source_states.get(unitig_id, "missing_dosage")
        if unitig_id in invalid:
            dosage, invalid_reason = invalid[unitig_id]
            candidate_decisions.append(
                AllelicCandidateDecision(
                    candidate,
                    source_state,
                    dosage,
                    None,
                    "unassigned",
                    invalid_reason,
                )
            )
            continue
        dosage = dosage_by_unitig[unitig_id]
        if dosage is None:
            raise AssertionError("eligible dosage unexpectedly missing")
        weighted_score = candidate.union_support_bases * dosage
        if status == "ambiguous" and unitig_id in optimal_members:
            candidate_status: CandidateStatus = "ambiguous"
            candidate_reason = "candidate_in_equal_optimum_subset"
        elif unitig_id in selected_unitigs:
            candidate_status = "selected"
            candidate_reason = "selected_unique_optimum"
        else:
            candidate_status = "rejected"
            candidate_reason = "lower_objective_under_capacity"
        candidate_decisions.append(
            AllelicCandidateDecision(
                candidate,
                source_state,
                dosage,
                weighted_score,
                candidate_status,
                candidate_reason,
            )
        )

    return AllelicBinDecision(
        bin_key=key,
        ploidy=ploidy,
        min_support_bases=min_support_bases,
        min_bin_coverage=min_bin_coverage,
        status=status,
        reason=reason,
        selected_unitigs=selected_unitigs,
        selected_copy_count=selected_copy_count,
        remaining_capacity=ploidy - selected_copy_count,
        best_objective=best_score,
        optimal_solution_count=optimal_count,
        candidates=tuple(candidate_decisions),
    )


def _merge_equal_states(
    first: _SolutionState,
    second: _SolutionState,
) -> _SolutionState:
    return _SolutionState(
        score=first.score,
        count=first.count + second.count,
        representative=min(first.representative, second.representative),
        member_union=first.member_union | second.member_union,
    )
