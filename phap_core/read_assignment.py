"""Deterministic, auditable assignment of reads to haplotype groups."""

from __future__ import annotations

import hashlib
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Dict, Optional, Tuple


@dataclass(frozen=True)
class UnitigCandidate:
    """Validated haplotype-group candidates for one unitig."""

    unitig_id: str
    dosage: Optional[int]
    groups: Tuple[str, ...]
    status: str
    reason: str


@dataclass(frozen=True)
class ReadAssignment:
    """One global assignment decision for a read or read pair."""

    read_id: str
    modality: str
    status: str
    destination_group: Optional[str]
    candidate_groups: Tuple[str, ...]
    unitigs: Tuple[str, ...]
    reason: str


def canonical_read_id(read_id: str, *, paired: bool) -> str:
    """Return one stable entity ID, merging common ``/1`` and ``/2`` mates."""
    if not read_id:
        raise ValueError("read_id must not be empty")
    if paired and read_id.endswith(("/1", "/2")):
        return read_id[:-2]
    return read_id


def paired_fastq_patterns(read_id: str, *, mate: int) -> Tuple[str, str]:
    """Return canonical and mate-suffixed patterns for common FASTQ IDs."""
    if mate not in {1, 2}:
        raise ValueError("mate must be 1 or 2")
    canonical = canonical_read_id(read_id, paired=True)
    return canonical, f"{canonical}/{mate}"


def alignment_is_usable(
    *,
    is_unmapped: bool,
    is_secondary: bool,
    is_supplementary: bool,
    is_duplicate: bool,
    is_qcfail: bool,
    is_proper_pair: bool,
    mapping_quality: int,
    min_mapq: int,
) -> bool:
    """Apply the read-assignment BAM policy.

    Unmapped, secondary, supplementary, duplicate, QC-fail, and low-MAPQ
    records are filtered. ``is_proper_pair`` is accepted for audit clarity but
    intentionally does not affect the result because Hi-C pairs are commonly
    non-proper under conventional paired-end alignment semantics.
    """
    if min_mapq < 0:
        raise ValueError("min_mapq must be non-negative")
    _ = is_proper_pair
    return not (
        is_unmapped
        or is_secondary
        or is_supplementary
        or is_duplicate
        or is_qcfail
        or mapping_quality < min_mapq
    )


def parse_unitig_dosages(
    lines: Iterable[str],
) -> Tuple[Dict[str, Optional[int]], Dict[str, str]]:
    """Parse current or legacy unitig-dosage tables without guessing invalid calls."""
    legacy_dosage = {
        "haplotig": 1,
        "diplotig": 2,
        "triplotig": 3,
        "tetraplotig": 4,
    }
    invalid_labels = {
        "ambiguous",
        "low_coverage",
        "high_copy",
        "mixed",
        "other",
    }
    rows = [line.strip().split() for line in lines if line.strip()]
    if not rows:
        raise ValueError("empty dosage table")

    header = rows[0]
    has_header = "contig_ID" in header
    data_rows = rows[1:] if has_header else rows
    id_index = header.index("contig_ID") if has_header else 0
    dosage_index = header.index("dosage") if has_header and "dosage" in header else None
    label_index: Optional[int] = None
    if has_header:
        for name in ("contig_type", "classification", "status", "dosage_class"):
            if name in header:
                label_index = header.index(name)
                break
    elif len(header) >= 3:
        label_index = 2

    dosage_by_unitig: Dict[str, Optional[int]] = {}
    source_labels: Dict[str, str] = {}
    for fields in data_rows:
        if id_index >= len(fields):
            raise ValueError(f"malformed dosage row: {' '.join(fields)}")
        unitig = fields[id_index]
        if not unitig:
            raise ValueError("empty unitig identifier in dosage table")
        if unitig in dosage_by_unitig:
            raise ValueError(f"duplicate unitig in dosage table: {unitig}")
        label = (
            fields[label_index].lower()
            if label_index is not None and label_index < len(fields)
            else ""
        )
        source_labels[unitig] = label or "numeric"

        dosage: Optional[int] = None
        if label not in invalid_labels and dosage_index is not None:
            if dosage_index < len(fields):
                try:
                    numeric = float(fields[dosage_index])
                    if numeric.is_integer():
                        dosage = int(numeric)
                except ValueError:
                    dosage = None
        elif label in legacy_dosage:
            dosage = legacy_dosage[label]
        elif label.startswith("dosage_"):
            try:
                dosage = int(label.removeprefix("dosage_"))
            except ValueError:
                dosage = None

        dosage_by_unitig[unitig] = dosage

    if not dosage_by_unitig:
        raise ValueError("dosage table contains no unitig records")
    return dosage_by_unitig, source_labels


def build_unitig_candidates(
    group_unitigs: Mapping[str, Iterable[str]],
    unitig_dosage: Mapping[str, Optional[int]],
    ploidy: int,
) -> Dict[str, UnitigCandidate]:
    """Validate dosage against group membership for every observed unitig."""
    if ploidy < 1:
        raise ValueError("ploidy must be at least 1")

    memberships: Dict[str, set[str]] = defaultdict(set)
    for group in sorted(group_unitigs):
        if not group:
            raise ValueError("group names must not be empty")
        for unitig in group_unitigs[group]:
            if unitig:
                memberships[unitig].add(group)

    candidates: Dict[str, UnitigCandidate] = {}
    all_unitigs = sorted(set(memberships) | set(unitig_dosage))
    for unitig in all_unitigs:
        groups = tuple(sorted(memberships.get(unitig, set())))
        dosage = unitig_dosage.get(unitig)
        if unitig not in unitig_dosage:
            status = "unassigned"
            reason = "missing_dosage"
        elif dosage is None:
            status = "unassigned"
            reason = "invalid_dosage_call"
        elif dosage < 1 or dosage > ploidy:
            status = "unassigned"
            reason = "dosage_outside_ploidy"
        elif not groups:
            status = "unassigned"
            reason = "no_candidate_group"
        elif len(groups) != dosage:
            status = "ambiguous"
            reason = "dosage_group_count_mismatch"
        else:
            status = "assignable"
            reason = "dosage_conserved"

        candidates[unitig] = UnitigCandidate(
            unitig_id=unitig,
            dosage=dosage,
            groups=groups,
            status=status,
            reason=reason,
        )
    return candidates


def stable_destination(
    *,
    read_id: str,
    modality: str,
    unitigs: Sequence[str],
    candidate_groups: Sequence[str],
    seed: int,
) -> str:
    """Choose one destination with a stable BLAKE2 hash."""
    groups = tuple(sorted(set(candidate_groups)))
    if not groups:
        raise ValueError("candidate_groups must not be empty")

    fields = [
        "PHap-read-assignment-v1",
        str(seed),
        modality,
        read_id,
        "\x1f".join(sorted(set(unitigs))),
    ]
    digest = hashlib.blake2b(
        "\x1e".join(fields).encode("utf-8"),
        digest_size=8,
        person=b"PHap.read.v1",
    ).digest()
    return groups[int.from_bytes(digest, byteorder="big") % len(groups)]


def assign_reads(
    read_unitigs: Mapping[str, Iterable[str]],
    unitig_candidates: Mapping[str, UnitigCandidate],
    *,
    modality: str,
    seed: int,
) -> Tuple[ReadAssignment, ...]:
    """Assign every read entity exactly once using compatible group evidence.

    Multi-mapping is resolved only when all mapped, assignable unitigs share at
    least one candidate group. Conflicting or invalid unitig evidence is kept
    out of haplotype FASTQ files and reported explicitly.
    """
    if not modality:
        raise ValueError("modality must not be empty")

    assignments: list[ReadAssignment] = []
    for read_id in sorted(read_unitigs):
        unitigs = tuple(sorted(set(read_unitigs[read_id])))
        if not unitigs:
            assignments.append(
                ReadAssignment(
                    read_id=read_id,
                    modality=modality,
                    status="unassigned",
                    destination_group=None,
                    candidate_groups=(),
                    unitigs=(),
                    reason="no_usable_alignment",
                )
            )
            continue

        observed = [unitig_candidates.get(unitig) for unitig in unitigs]
        assignable = [
            candidate
            for candidate in observed
            if candidate is not None and candidate.status == "assignable"
        ]
        invalid = [
            candidate
            for candidate in observed
            if candidate is None or candidate.status != "assignable"
        ]
        all_groups = tuple(
            sorted(
                {
                    group
                    for candidate in observed
                    if candidate is not None
                    for group in candidate.groups
                }
            )
        )

        if invalid:
            invalid_statuses = {
                "unassigned" if candidate is None else candidate.status
                for candidate in invalid
            }
            status = "ambiguous" if assignable or "ambiguous" in invalid_statuses else "unassigned"
            reason = (
                "mixed_assignable_and_invalid_unitigs"
                if assignable
                else "invalid_unitig_candidates"
            )
            assignments.append(
                ReadAssignment(
                    read_id=read_id,
                    modality=modality,
                    status=status,
                    destination_group=None,
                    candidate_groups=all_groups,
                    unitigs=unitigs,
                    reason=reason,
                )
            )
            continue

        compatible_groups = set(assignable[0].groups)
        for candidate in assignable[1:]:
            compatible_groups.intersection_update(candidate.groups)
        compatible = tuple(sorted(compatible_groups))

        if not compatible:
            assignments.append(
                ReadAssignment(
                    read_id=read_id,
                    modality=modality,
                    status="ambiguous",
                    destination_group=None,
                    candidate_groups=all_groups,
                    unitigs=unitigs,
                    reason="incompatible_multimap",
                )
            )
            continue

        destination = stable_destination(
            read_id=read_id,
            modality=modality,
            unitigs=unitigs,
            candidate_groups=compatible,
            seed=seed,
        )
        if len(unitigs) > 1:
            reason = (
                "compatible_multimap_resolved"
                if len(compatible) == 1
                else "compatible_multimap_partition"
            )
        else:
            reason = "single_group" if len(compatible) == 1 else "collapsed_unitig_partition"
        assignments.append(
            ReadAssignment(
                read_id=read_id,
                modality=modality,
                status="assigned",
                destination_group=destination,
                candidate_groups=compatible,
                unitigs=unitigs,
                reason=reason,
            )
        )

    result = tuple(assignments)
    validate_assignments(result)
    return result


def validate_assignments(assignments: Sequence[ReadAssignment]) -> None:
    """Validate completeness and mutual exclusion of assignment decisions."""
    seen: set[Tuple[str, str]] = set()
    for assignment in assignments:
        entity = (assignment.modality, assignment.read_id)
        if entity in seen:
            raise ValueError(
                f"duplicate assignment decision for {assignment.modality}:{assignment.read_id}"
            )
        seen.add(entity)
        if assignment.status == "assigned":
            if assignment.destination_group is None:
                raise ValueError("assigned reads must have a destination group")
            if assignment.destination_group not in assignment.candidate_groups:
                raise ValueError("destination group must be one of the candidate groups")
        elif assignment.destination_group is not None:
            raise ValueError("ambiguous/unassigned reads must not have a destination group")


def group_assigned_reads(
    assignments: Sequence[ReadAssignment],
    groups: Iterable[str],
) -> Dict[str, Tuple[str, ...]]:
    """Return sorted disjoint read IDs for every declared group."""
    validate_assignments(assignments)
    grouped: Dict[str, list[str]] = {
        group: [] for group in sorted(set(groups))
    }
    for assignment in assignments:
        if assignment.status != "assigned":
            continue
        destination = assignment.destination_group
        if destination is None or destination not in grouped:
            raise ValueError(f"unknown destination group: {destination}")
        grouped[destination].append(assignment.read_id)

    read_owner: Dict[str, str] = {}
    result: Dict[str, Tuple[str, ...]] = {}
    for group in sorted(grouped):
        read_ids = tuple(sorted(grouped[group]))
        for read_id in read_ids:
            if read_id in read_owner:
                raise ValueError(
                    f"read {read_id} assigned to both {read_owner[read_id]} and {group}"
                )
            read_owner[read_id] = group
        result[group] = read_ids
    return result
