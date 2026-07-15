"""Atomic preparation of chromosome-specific clustering tables."""

from __future__ import annotations

from pathlib import Path
from typing import Union

from .atomic_io import atomic_write_lines


def write_chromosome_allelic_table(
    source_path: Union[str, Path],
    chromosome: str,
    destination_path: Union[str, Path],
) -> int:
    """Atomically replace one chromosome table from the current global table."""

    if not chromosome or "\t" in chromosome or "\n" in chromosome:
        raise ValueError("chromosome must be a non-empty tab-free identifier")
    source = Path(source_path)
    destination = Path(destination_path)
    if source.resolve() == destination.resolve():
        raise ValueError("source and destination allelic tables must differ")

    selected: list[str] = []
    with source.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.rstrip("\r\n")
            if not text:
                continue
            fields = text.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"{source}:{line_number}: expected chromosome, start, end, unitigs"
                )
            if fields[0] == chromosome:
                selected.append(text)
    atomic_write_lines(destination, selected)
    return len(selected)
