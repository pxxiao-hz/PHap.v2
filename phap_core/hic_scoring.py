"""Deterministic Hi-C scoring for dosage-aware haplotype-group assignment.

The normalized score for candidate unitig ``u`` and group ``g`` is

``sum(links(u, v) for v in g) / sum(RE_sites(v) for v in g)``.

The candidate's own restriction-site count is intentionally omitted. It is a
constant when ranking groups for the same candidate and therefore cannot
change that ranking. As a consequence, ``re_density`` values are only
comparable between groups for one candidate, not between different candidate
unitigs.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from numbers import Real
from typing import Literal, Optional, Tuple


ScoreMode = Literal["raw", "re_density"]


@dataclass(frozen=True)
class HiCGroupScore:
    """Auditable raw and restriction-site-normalized evidence for one group."""

    group_id: str
    raw_score: float
    re_sites: Optional[int]
    re_density: Optional[float]
    re_status: str
    raw_rank: int
    re_density_rank: Optional[int]


@dataclass(frozen=True)
class HiCGroupAssignment:
    """One dosage-aware group-assignment decision for a candidate unitig."""

    unitig_id: str
    dosage: int
    ploidy: int
    mode: ScoreMode
    status: str
    selected_groups: Tuple[str, ...]
    reason: str
    scores: Tuple[HiCGroupScore, ...]
    top_margin: Optional[float]
    selection_margin: Optional[float]
    min_score: float
    min_margin: float


@dataclass(frozen=True)
class _UnrankedScore:
    group_id: str
    raw_score: float
    re_sites: Optional[int]
    re_density: Optional[float]
    re_status: str


def validate_pair_links(
    pair_links: Mapping[Tuple[str, str], float],
) -> None:
    """Validate one-count-per-unordered-pair Hi-C input.

    A mapping cannot contain the same oriented key twice, but it can contain
    both ``(u, v)`` and ``(v, u)``. Accepting both would silently count the
    same undirected contig pair twice during candidate scoring.
    """

    seen_unordered: set[Tuple[str, str]] = set()
    for pair, count in pair_links.items():
        if (
            not isinstance(pair, tuple)
            or len(pair) != 2
            or not isinstance(pair[0], str)
            or not isinstance(pair[1], str)
            or not pair[0]
            or not pair[1]
        ):
            raise ValueError(
                "Hi-C pair keys must be two-element tuples of non-empty unitig IDs"
            )
        if isinstance(count, bool) or not isinstance(count, Real):
            raise ValueError("Hi-C link counts must be real numbers")
        if not math.isfinite(count) or count < 0:
            raise ValueError("Hi-C link counts must be finite and non-negative")
        unordered = tuple(sorted(pair))
        if unordered in seen_unordered:
            raise ValueError(
                f"Hi-C pair {unordered!r} occurs in both orientations"
            )
        seen_unordered.add(unordered)


def calculate_group_scores(
    unitig_id: str,
    group_unitigs: Mapping[str, Iterable[str]],
    links_by_unitig: Mapping[str, float],
    restriction_sites: Mapping[str, int],
) -> Tuple[HiCGroupScore, ...]:
    """Calculate raw and normalized scores for every candidate group.

    ``links_by_unitig[v]`` must contain the already filtered count of Hi-C
    pairs between the candidate and member unitig ``v``. A missing entry means
    zero observed links. Pair filtering and pair de-duplication belong at the
    BAM parsing boundary, before this pure scoring function.

    Restriction-site counts must be available for every member of a group.
    Missing counts and a group total of zero are explicitly unsupported; the
    denominator is never replaced with one.
    """
    if not unitig_id:
        raise ValueError("unitig_id must not be empty")
    if not group_unitigs:
        raise ValueError("group_unitigs must not be empty")

    _validate_link_counts(links_by_unitig)
    _validate_restriction_sites(restriction_sites)

    unranked = []
    for group_id in sorted(group_unitigs):
        if not group_id:
            raise ValueError("group IDs must not be empty")
        members = tuple(sorted(set(group_unitigs[group_id])))
        if any(not member for member in members):
            raise ValueError(f"group {group_id!r} contains an empty unitig ID")
        if unitig_id in members:
            raise ValueError(
                f"candidate unitig {unitig_id!r} must not already belong to group {group_id!r}"
            )

        raw_score = math.fsum(links_by_unitig.get(member, 0.0) for member in members)
        missing_re = tuple(member for member in members if member not in restriction_sites)
        if missing_re:
            re_sites: Optional[int] = None
            re_density: Optional[float] = None
            re_status = "missing_re_sites"
        else:
            re_sites = sum(restriction_sites[member] for member in members)
            if re_sites == 0:
                re_density = None
                re_status = "zero_re_sites"
            else:
                re_density = raw_score / re_sites
                re_status = "supported"

        unranked.append(
            _UnrankedScore(
                group_id=group_id,
                raw_score=raw_score,
                re_sites=re_sites,
                re_density=re_density,
                re_status=re_status,
            )
        )

    raw_ranks = {
        score.group_id: rank
        for rank, score in enumerate(
            sorted(unranked, key=lambda item: (-item.raw_score, item.group_id)),
            start=1,
        )
    }
    density_ranks = {
        score.group_id: rank
        for rank, score in enumerate(
            sorted(
                (item for item in unranked if item.re_density is not None),
                key=lambda item: (-_require_density(item), item.group_id),
            ),
            start=1,
        )
    }
    return tuple(
        HiCGroupScore(
            group_id=score.group_id,
            raw_score=score.raw_score,
            re_sites=score.re_sites,
            re_density=score.re_density,
            re_status=score.re_status,
            raw_rank=raw_ranks[score.group_id],
            re_density_rank=density_ranks.get(score.group_id),
        )
        for score in unranked
    )


def assign_candidate_groups(
    unitig_id: str,
    group_unitigs: Mapping[str, Iterable[str]],
    links_by_unitig: Mapping[str, float],
    restriction_sites: Mapping[str, int],
    *,
    dosage: int,
    ploidy: int,
    mode: ScoreMode,
    min_score: float,
    min_margin: float,
) -> HiCGroupAssignment:
    """Select exactly ``dosage`` groups or return an explicit ambiguity.

    ``min_score`` is measured in raw pair counts for ``mode="raw"`` and in
    links per group restriction site for ``mode="re_density"``.
    ``min_margin`` applies to the difference between the last selected and
    first rejected group. Exact boundary ties are always ambiguous, including
    when ``min_margin`` is zero. Ties wholly inside the selected set are safe
    because they do not change membership.

    When ``dosage == ploidy``, every group is selected by copy-count
    conservation and no ranking boundary exists. Scores remain in the audit
    object, including unsupported RE normalization.
    """
    _validate_assignment_parameters(
        group_unitigs=group_unitigs,
        dosage=dosage,
        ploidy=ploidy,
        mode=mode,
        min_score=min_score,
        min_margin=min_margin,
    )
    scores = calculate_group_scores(
        unitig_id,
        group_unitigs,
        links_by_unitig,
        restriction_sites,
    )
    ranked = _rank_for_mode(scores, mode)
    top_margin = _margin(ranked, 1)
    selection_margin = _margin(ranked, dosage)

    if dosage == ploidy:
        return HiCGroupAssignment(
            unitig_id=unitig_id,
            dosage=dosage,
            ploidy=ploidy,
            mode=mode,
            status="assigned",
            selected_groups=tuple(sorted(group_unitigs)),
            reason="dosage_equals_ploidy",
            scores=scores,
            top_margin=top_margin,
            selection_margin=None,
            min_score=min_score,
            min_margin=min_margin,
        )

    if len(ranked) != ploidy:
        return _ambiguous_assignment(
            unitig_id=unitig_id,
            dosage=dosage,
            ploidy=ploidy,
            mode=mode,
            reason="unsupported_re_density",
            scores=scores,
            top_margin=top_margin,
            selection_margin=selection_margin,
            min_score=min_score,
            min_margin=min_margin,
        )

    weakest_selected_score = ranked[dosage - 1][1]
    if weakest_selected_score < min_score:
        return _ambiguous_assignment(
            unitig_id=unitig_id,
            dosage=dosage,
            ploidy=ploidy,
            mode=mode,
            reason="below_min_score",
            scores=scores,
            top_margin=top_margin,
            selection_margin=selection_margin,
            min_score=min_score,
            min_margin=min_margin,
        )

    if selection_margin is None:
        raise AssertionError("a partial-ploidy assignment must have a ranking boundary")
    if selection_margin == 0.0:
        reason = "non_unique_boundary"
    elif selection_margin < min_margin:
        reason = "below_min_margin"
    else:
        reason = ""
    if reason:
        return _ambiguous_assignment(
            unitig_id=unitig_id,
            dosage=dosage,
            ploidy=ploidy,
            mode=mode,
            reason=reason,
            scores=scores,
            top_margin=top_margin,
            selection_margin=selection_margin,
            min_score=min_score,
            min_margin=min_margin,
        )

    return HiCGroupAssignment(
        unitig_id=unitig_id,
        dosage=dosage,
        ploidy=ploidy,
        mode=mode,
        status="assigned",
        selected_groups=tuple(sorted(group_id for group_id, _ in ranked[:dosage])),
        reason="dosage_conserved",
        scores=scores,
        top_margin=top_margin,
        selection_margin=selection_margin,
        min_score=min_score,
        min_margin=min_margin,
    )


def _validate_link_counts(links_by_unitig: Mapping[str, float]) -> None:
    for unitig_id, count in links_by_unitig.items():
        if not unitig_id:
            raise ValueError("links_by_unitig contains an empty unitig ID")
        if not math.isfinite(count) or count < 0:
            raise ValueError(f"Hi-C link count for {unitig_id!r} must be finite and non-negative")


def _validate_restriction_sites(restriction_sites: Mapping[str, int]) -> None:
    for unitig_id, count in restriction_sites.items():
        if not unitig_id:
            raise ValueError("restriction_sites contains an empty unitig ID")
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ValueError(
                f"restriction-site count for {unitig_id!r} must be a non-negative integer"
            )


def _validate_assignment_parameters(
    *,
    group_unitigs: Mapping[str, Iterable[str]],
    dosage: int,
    ploidy: int,
    mode: ScoreMode,
    min_score: float,
    min_margin: float,
) -> None:
    if isinstance(ploidy, bool) or not isinstance(ploidy, int) or ploidy < 1:
        raise ValueError("ploidy must be a positive integer")
    if isinstance(dosage, bool) or not isinstance(dosage, int):
        raise ValueError("dosage must be an integer")
    if not 1 <= dosage <= ploidy:
        raise ValueError("dosage must be between 1 and ploidy")
    if len(group_unitigs) != ploidy:
        raise ValueError("group_unitigs must contain exactly ploidy groups")
    if mode not in ("raw", "re_density"):
        raise ValueError("mode must be 'raw' or 're_density'")
    if not math.isfinite(min_score) or min_score < 0:
        raise ValueError("min_score must be finite and non-negative")
    if not math.isfinite(min_margin) or min_margin < 0:
        raise ValueError("min_margin must be finite and non-negative")


def _require_density(score: _UnrankedScore) -> float:
    if score.re_density is None:
        raise AssertionError("density ranking received an unsupported score")
    return score.re_density


def _rank_for_mode(
    scores: Tuple[HiCGroupScore, ...], mode: ScoreMode
) -> Tuple[Tuple[str, float], ...]:
    if mode == "raw":
        values = ((score.group_id, score.raw_score) for score in scores)
    else:
        values = (
            (score.group_id, score.re_density)
            for score in scores
            if score.re_density is not None
        )
    return tuple(sorted(values, key=lambda item: (-item[1], item[0])))


def _margin(ranked: Tuple[Tuple[str, float], ...], boundary: int) -> Optional[float]:
    if boundary < 1:
        raise ValueError("boundary must be positive")
    if len(ranked) <= boundary:
        return None
    return ranked[boundary - 1][1] - ranked[boundary][1]


def _ambiguous_assignment(
    *,
    unitig_id: str,
    dosage: int,
    ploidy: int,
    mode: ScoreMode,
    reason: str,
    scores: Tuple[HiCGroupScore, ...],
    top_margin: Optional[float],
    selection_margin: Optional[float],
    min_score: float,
    min_margin: float,
) -> HiCGroupAssignment:
    return HiCGroupAssignment(
        unitig_id=unitig_id,
        dosage=dosage,
        ploidy=ploidy,
        mode=mode,
        status="ambiguous",
        selected_groups=(),
        reason=reason,
        scores=scores,
        top_margin=top_margin,
        selection_margin=selection_margin,
        min_score=min_score,
        min_margin=min_margin,
    )
