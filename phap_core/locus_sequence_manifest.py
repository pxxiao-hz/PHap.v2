"""Validation of the current locus FASTA publication manifest."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Union

from .fasta import count_fasta_records


PathLike = Union[str, Path]
LOCUS_FASTA_SUFFIX = ".putg.fa"
_SAFE_LOCUS_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


@dataclass(frozen=True)
class LocusSequenceOutput:
    """One current-locus FASTA declared by the step-2 manifest."""

    locus_id: str
    path: Path
    record_count: int


def read_current_locus_outputs(
    output_directory: PathLike,
) -> tuple[LocusSequenceOutput, ...]:
    """Read and fully validate the current step-2 locus FASTA publication."""

    root = Path(output_directory).resolve()
    manifest = root / "locus_sequence_manifest.tsv"
    if not manifest.is_file():
        raise ValueError(f"locus sequence manifest is missing: {manifest}")
    rows = manifest.read_text(encoding="utf-8").splitlines()
    if not rows or rows[0] != "output_file\trecord_count":
        raise ValueError(f"invalid locus sequence manifest header: {manifest}")

    declared_names: set[str] = set()
    outputs: list[LocusSequenceOutput] = []
    unplaced_seen = False
    for line_number, line in enumerate(rows[1:], start=2):
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError(f"{manifest}:{line_number}: expected two columns")
        filename, count_text = fields
        if not filename or Path(filename).name != filename:
            raise ValueError(f"{manifest}:{line_number}: unsafe output filename")
        if filename in declared_names:
            raise ValueError(f"{manifest}:{line_number}: duplicate output filename")
        declared_names.add(filename)
        try:
            record_count = int(count_text)
        except ValueError as error:
            raise ValueError(
                f"{manifest}:{line_number}: invalid FASTA record count"
            ) from error
        if record_count < 0:
            raise ValueError(f"{manifest}:{line_number}: negative FASTA record count")

        path = (root / filename).resolve()
        try:
            path.relative_to(root)
        except ValueError as error:
            raise ValueError(f"{manifest}:{line_number}: output escapes its root") from error
        if filename == "un_chr.fa":
            unplaced_seen = True
            _validate_fasta_count(path, record_count)
            continue

        locus_id = locus_id_from_fasta_name(filename)
        if record_count == 0:
            raise ValueError(f"{manifest}:{line_number}: current locus FASTA is empty")
        _validate_fasta_count(path, record_count)
        outputs.append(LocusSequenceOutput(locus_id, path, record_count))

    if not unplaced_seen:
        raise ValueError(f"locus sequence manifest omits un_chr.fa: {manifest}")
    existing_locus_files = {path.name for path in root.glob(f"*{LOCUS_FASTA_SUFFIX}")}
    declared_locus_files = {
        output.path.name
        for output in outputs
    }
    if existing_locus_files != declared_locus_files:
        unexpected = sorted(existing_locus_files - declared_locus_files)
        missing = sorted(declared_locus_files - existing_locus_files)
        details = []
        if unexpected:
            details.append("unexpected=" + ",".join(unexpected))
        if missing:
            details.append("missing=" + ",".join(missing))
        raise ValueError("locus FASTA set does not match manifest: " + "; ".join(details))

    ordered = tuple(sorted(outputs, key=lambda row: row.locus_id))
    if len({row.locus_id for row in ordered}) != len(ordered):
        raise ValueError("locus sequence manifest contains duplicate locus IDs")
    return ordered


def locus_id_from_fasta_name(filename: str) -> str:
    """Remove one terminal locus FASTA suffix without altering the locus ID."""

    if Path(filename).name != filename or not filename.endswith(LOCUS_FASTA_SUFFIX):
        raise ValueError(f"invalid locus FASTA filename: {filename!r}")
    locus_id = filename[: -len(LOCUS_FASTA_SUFFIX)]
    if _SAFE_LOCUS_ID.fullmatch(locus_id) is None:
        raise ValueError(f"unsafe or empty locus ID: {locus_id!r}")
    return locus_id


def _validate_fasta_count(path: Path, expected_count: int) -> None:
    if not path.is_file():
        raise ValueError(f"manifest-declared FASTA is missing: {path}")
    if expected_count == 0:
        if path.stat().st_size != 0:
            raise ValueError(f"FASTA record count mismatch for {path}: expected 0")
        return
    observed_count = count_fasta_records(path)
    if observed_count != expected_count:
        raise ValueError(
            f"FASTA record count mismatch for {path}: expected {expected_count}, "
            f"observed {observed_count}"
        )
