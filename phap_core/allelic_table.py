"""Dosage-safe filtering of mT2T allelic-table rows."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass(frozen=True)
class AllelicTableRow:
    chromosome: str
    start: int
    end: int
    unitigs: Tuple[str, ...]


@dataclass(frozen=True)
class AllelicExclusion:
    unitig_id: str
    chromosome: str
    start: int
    end: int
    source_state: str
    reason: str


def find_missing_bridge_unitigs(
    rows: Sequence[AllelicTableRow],
    index: int,
    *,
    search_range: int,
) -> Tuple[str, ...]:
    """Find same-chromosome unitigs spanning a gap in one allelic bin."""

    if not 0 <= index < len(rows):
        raise IndexError("allelic-table row index is out of range")
    if search_range < 0:
        raise ValueError("search_range must be non-negative")
    current = rows[index]
    _validate_row(current)
    previous: set[str] = set()
    for row in rows[max(0, index - search_range) : index]:
        if row.chromosome == current.chromosome:
            previous.update(row.unitigs)
    following: set[str] = set()
    for row in rows[index + 1 : index + 1 + search_range]:
        if row.chromosome == current.chromosome:
            following.update(row.unitigs)
    return tuple(sorted((previous & following) - set(current.unitigs)))


def filter_allelic_rows(
    rows: Sequence[AllelicTableRow],
    dosage_by_unitig: Mapping[str, Optional[int]],
    source_states: Mapping[str, str],
    unitig_lengths: Mapping[str, int],
    *,
    ploidy: int,
) -> Tuple[Tuple[Optional[AllelicTableRow], ...], Tuple[AllelicExclusion, ...]]:
    """Remove invalid dosage calls and reject unresolved over-capacity bins."""

    if ploidy < 1:
        raise ValueError("ploidy must be positive")
    seen_unitigs: set[str] = set()
    filtered: list[Optional[AllelicTableRow]] = []
    audit: list[AllelicExclusion] = []
    total_bins = len(rows)

    for index, row in enumerate(rows):
        _validate_row(row)
        ordered = sorted(
            set(row.unitigs),
            key=lambda unitig_id: (-unitig_lengths.get(unitig_id, 0), unitig_id),
        )
        retained = []
        repeated = []
        copy_count = 0
        for unitig_id in ordered:
            dosage = dosage_by_unitig.get(unitig_id)
            if dosage is None or not 1 <= dosage <= ploidy:
                audit.append(
                    _exclusion(
                        row,
                        unitig_id,
                        source_states,
                        "invalid_or_missing_dosage",
                    )
                )
                continue
            retained.append(unitig_id)
            copy_count += dosage
            if unitig_id in seen_unitigs:
                repeated.append(unitig_id)

        if copy_count > ploidy and 0 < index < total_bins - 1:
            previous = set(rows[index - 1].unitigs)
            following = set(rows[index + 1].unitigs)
            new_unitigs = set(row.unitigs) - previous
            if new_unitigs and all(unitig in following for unitig in new_unitigs):
                for unitig_id in retained:
                    audit.append(
                        _exclusion(
                            row,
                            unitig_id,
                            source_states,
                            "redundant_over_capacity_bin",
                        )
                    )
                filtered.append(None)
                continue

        while copy_count > ploidy and repeated:
            unitig_id = repeated.pop()
            retained.remove(unitig_id)
            dosage = dosage_by_unitig.get(unitig_id)
            if dosage is None:
                raise AssertionError("validated dosage unexpectedly missing")
            copy_count -= dosage
            audit.append(
                _exclusion(
                    row,
                    unitig_id,
                    source_states,
                    "repeated_unitig_exceeds_bin_capacity",
                )
            )
        if copy_count > ploidy:
            for unitig_id in retained:
                audit.append(
                    _exclusion(
                        row,
                        unitig_id,
                        source_states,
                        "ambiguous_bin_exceeds_ploidy",
                    )
                )
            filtered.append(None)
            continue

        retained_tuple = tuple(sorted(retained))
        if not retained_tuple:
            filtered.append(None)
            continue
        filtered.append(
            AllelicTableRow(row.chromosome, row.start, row.end, retained_tuple)
        )
        seen_unitigs.update(retained_tuple)

    return tuple(filtered), tuple(
        sorted(
            audit,
            key=lambda row: (
                row.chromosome,
                row.start,
                row.end,
                row.unitig_id,
                row.reason,
            ),
        )
    )


def _validate_row(row: AllelicTableRow) -> None:
    if not row.chromosome:
        raise ValueError("chromosome must not be empty")
    if row.start < 0 or row.end < row.start:
        raise ValueError("allelic-table coordinates must be 0-based half-open")
    if not row.unitigs or any(not unitig_id for unitig_id in row.unitigs):
        raise ValueError("allelic-table rows must contain non-empty unitig IDs")


def _exclusion(
    row: AllelicTableRow,
    unitig_id: str,
    source_states: Mapping[str, str],
    reason: str,
) -> AllelicExclusion:
    return AllelicExclusion(
        unitig_id,
        row.chromosome,
        row.start,
        row.end,
        source_states.get(unitig_id, "missing_dosage"),
        reason,
    )
