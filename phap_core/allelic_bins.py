"""Pure construction and ranking of mT2T allelic-bin candidates."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import DefaultDict, Dict, Tuple

from .locus_rescue import interval_union_length


@dataclass(frozen=True, order=True)
class AllelicBinKey:
    """One 0-based half-open interval on an mT2T target."""

    target_id: str
    start: int
    end: int


@dataclass(frozen=True)
class TargetIntervalEvidence:
    """One validated target interval contributed by a unitig alignment."""

    unitig_id: str
    target_id: str
    target_length: int
    target_start: int
    target_end: int


@dataclass(frozen=True)
class AllelicBinCandidate:
    """Interval-union support for one unitig in one target bin."""

    bin_key: AllelicBinKey
    unitig_id: str
    union_support_bases: int


def build_allelic_bin_candidates(
    alignments: Iterable[TargetIntervalEvidence],
    *,
    bin_size: int,
) -> Tuple[AllelicBinCandidate, ...]:
    """Split target intervals into bins and union overlaps per unitig."""

    if bin_size < 1:
        raise ValueError("bin_size must be positive")

    target_lengths: Dict[str, int] = {}
    intervals: DefaultDict[
        Tuple[AllelicBinKey, str], list[Tuple[int, int]]
    ] = defaultdict(list)
    for alignment in alignments:
        _validate_alignment(alignment)
        previous_length = target_lengths.setdefault(
            alignment.target_id,
            alignment.target_length,
        )
        if previous_length != alignment.target_length:
            raise ValueError(
                f"conflicting lengths for target {alignment.target_id!r}: "
                f"{previous_length} and {alignment.target_length}"
            )

        position = alignment.target_start
        while position < alignment.target_end:
            bin_start = (position // bin_size) * bin_size
            bin_end = min(bin_start + bin_size, alignment.target_length)
            segment_end = min(alignment.target_end, bin_end)
            key = AllelicBinKey(alignment.target_id, bin_start, bin_end)
            intervals[(key, alignment.unitig_id)].append((position, segment_end))
            position = segment_end

    return tuple(
        AllelicBinCandidate(
            bin_key=key,
            unitig_id=unitig_id,
            union_support_bases=interval_union_length(
                candidate_intervals,
                sequence_length=key.end,
            ),
        )
        for (key, unitig_id), candidate_intervals in sorted(intervals.items())
    )


def select_legacy_top_n_candidates(
    candidates: Sequence[AllelicBinCandidate],
    *,
    top_n: int,
) -> Tuple[AllelicBinCandidate, ...]:
    """Reproduce the historical count-based per-bin ranking deterministically."""

    if top_n < 1:
        raise ValueError("top_n must be positive")

    by_bin: DefaultDict[AllelicBinKey, list[AllelicBinCandidate]] = defaultdict(list)
    seen: set[Tuple[AllelicBinKey, str]] = set()
    for candidate in candidates:
        _validate_candidate(candidate)
        identity = (candidate.bin_key, candidate.unitig_id)
        if identity in seen:
            raise ValueError(
                "duplicate allelic-bin candidate for "
                f"{candidate.bin_key.target_id}:{candidate.bin_key.start}-"
                f"{candidate.bin_key.end}:{candidate.unitig_id}"
            )
        seen.add(identity)
        by_bin[candidate.bin_key].append(candidate)

    selected: list[AllelicBinCandidate] = []
    for key in sorted(by_bin):
        ranked = sorted(
            by_bin[key],
            key=lambda candidate: (
                -candidate.union_support_bases,
                candidate.unitig_id,
            ),
        )
        selected.extend(ranked[:top_n])
    return tuple(selected)


def _validate_alignment(alignment: TargetIntervalEvidence) -> None:
    if not alignment.unitig_id or not alignment.target_id:
        raise ValueError("unitig_id and target_id must not be empty")
    if alignment.target_length < 1:
        raise ValueError("target_length must be positive")
    if not 0 <= alignment.target_start < alignment.target_end <= alignment.target_length:
        raise ValueError(
            "target interval must be non-empty, 0-based half-open, and within "
            "the declared target length"
        )


def _validate_candidate(candidate: AllelicBinCandidate) -> None:
    key = candidate.bin_key
    if not key.target_id or key.start < 0 or key.end <= key.start:
        raise ValueError("invalid allelic bin key")
    if not candidate.unitig_id:
        raise ValueError("unitig_id must not be empty")
    if not 0 < candidate.union_support_bases <= key.end - key.start:
        raise ValueError("union support must be within the bin length")
