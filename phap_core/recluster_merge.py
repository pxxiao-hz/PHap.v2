"""Deterministic merging of current-locus reclustering outputs."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, Union

from .atomic_io import atomic_text_writer, atomic_write_lines
from .stage_manifest import sha256_file


PathLike = Union[str, Path]
_SAFE_LOCUS_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


@dataclass(frozen=True)
class ReclusterMergeSource:
    """One verified current-locus cluster file included in the merge."""

    locus_id: str
    path: Path
    size: int
    sha256: str
    group_ids: tuple[str, ...]


def merge_current_recluster_outputs(
    recluster_root: PathLike,
    locus_ids: Sequence[str],
    destination_path: PathLike,
    audit_path: PathLike,
) -> tuple[ReclusterMergeSource, ...]:
    """Atomically merge only outputs belonging to the current locus set."""

    root = Path(recluster_root).resolve()
    destination = Path(destination_path).resolve()
    audit = Path(audit_path).resolve()
    if destination == audit:
        raise ValueError("merged cluster output and audit paths must differ")

    ordered_loci = tuple(sorted(locus_ids))
    if len(ordered_loci) != len(set(ordered_loci)):
        raise ValueError("current locus IDs must be unique")
    for locus_id in ordered_loci:
        if _SAFE_LOCUS_ID.fullmatch(locus_id) is None:
            raise ValueError(f"unsafe current locus ID: {locus_id!r}")

    sources = tuple(_capture_source(root, locus_id) for locus_id in ordered_loci)
    group_ids = tuple(group_id for source in sources for group_id in source.group_ids)
    if len(group_ids) != len(set(group_ids)):
        raise ValueError("current-locus recluster outputs contain duplicate groups")
    protected_paths = {root, *(source.path for source in sources)}
    if destination in protected_paths or audit in protected_paths:
        raise ValueError("merge outputs must not replace the recluster root or inputs")

    audit.unlink(missing_ok=True)
    with atomic_text_writer(destination) as merged:
        for source in sources:
            with source.path.open(encoding="utf-8", newline="") as handle:
                for line in handle:
                    merged.write(line)
        _require_unchanged_sources(sources)

    atomic_write_lines(
        audit,
        (
            "locus_ID\tsource_file\tsize\tsha256\tdestination\treason",
            *(
                "\t".join(
                    (
                        source.locus_id,
                        str(source.path.relative_to(root)),
                        str(source.size),
                        source.sha256,
                        destination.name,
                        "current_locus_recluster_output",
                    )
                )
                for source in sources
            ),
        ),
    )
    return sources


def _capture_source(root: Path, locus_id: str) -> ReclusterMergeSource:
    source = (root / locus_id / "group.reassignment.cluster.txt").resolve()
    try:
        source.relative_to(root)
    except ValueError as error:
        raise ValueError(f"recluster input escapes its root: {source}") from error
    if not source.is_file():
        raise ValueError(f"current-locus recluster output is missing: {source}")
    size = source.stat().st_size
    if size and not _ends_with_newline(source):
        raise ValueError(f"recluster output is not newline-terminated: {source}")
    group_ids = _parse_group_ids(source, locus_id)
    return ReclusterMergeSource(
        locus_id,
        source,
        size,
        sha256_file(source),
        group_ids,
    )


def _ends_with_newline(path: Path) -> bool:
    with path.open("rb") as handle:
        handle.seek(-1, 2)
        return handle.read(1) == b"\n"


def _parse_group_ids(path: Path, locus_id: str) -> tuple[str, ...]:
    group_ids: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            group_id = fields[0]
            if not group_id.startswith(f"{locus_id}_"):
                raise ValueError(
                    f"{path}:{line_number}: group {group_id!r} does not belong "
                    f"to locus {locus_id!r}"
                )
            if group_id in group_ids:
                raise ValueError(f"{path}:{line_number}: duplicate group {group_id!r}")
            unitigs = fields[1:]
            if len(unitigs) != len(set(unitigs)):
                raise ValueError(
                    f"{path}:{line_number}: duplicate unitig within group {group_id!r}"
                )
            group_ids.append(group_id)
    if not group_ids:
        raise ValueError(f"current-locus recluster output contains no groups: {path}")
    return tuple(group_ids)


def _require_unchanged_sources(sources: Sequence[ReclusterMergeSource]) -> None:
    for source in sources:
        if source.path.stat().st_size != source.size:
            raise ValueError(f"recluster input changed during merge: {source.path}")
        if sha256_file(source.path) != source.sha256:
            raise ValueError(f"recluster input changed during merge: {source.path}")
