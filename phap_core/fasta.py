"""Validated FASTA records with exact header preservation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Tuple, Union


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA record.

    ``identifier`` is the first whitespace-delimited header token used by PAF
    and BAM. ``header`` preserves the complete source header without ``>``.
    """

    identifier: str
    header: str
    sequence: str


def read_fasta(path: Union[str, Path]) -> Tuple[FastaRecord, ...]:
    """Read FASTA records, rejecting malformed or duplicate identifiers."""

    records: list[FastaRecord] = []
    identifiers: set[str] = set()
    current_header: Optional[str] = None
    sequence_parts: list[str] = []
    last_line_number = 0
    source = str(path)
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            last_line_number = line_number
            text = line.rstrip("\r\n")
            if not text:
                continue
            if text.startswith(">"):
                if current_header is not None:
                    records.append(
                        _make_record(
                            current_header,
                            sequence_parts,
                            identifiers,
                            source,
                            line_number - 1,
                        )
                    )
                current_header = text[1:]
                sequence_parts = []
            elif current_header is None:
                raise ValueError(
                    f"{source}:{line_number}: sequence appears before a FASTA header"
                )
            else:
                sequence_parts.append(text.strip())
    if current_header is not None:
        records.append(
            _make_record(
                current_header,
                sequence_parts,
                identifiers,
                source,
                last_line_number,
            )
        )
    if not records:
        raise ValueError(f"{source}: FASTA contains no records")
    return tuple(records)


def count_fasta_records(path: Union[str, Path]) -> int:
    """Validate a FASTA and count records without materializing sequences."""

    source = str(path)
    identifiers: set[str] = set()
    current_id: Optional[str] = None
    current_has_sequence = False
    count = 0
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.rstrip("\r\n")
            if not text:
                continue
            if text.startswith(">"):
                if current_id is not None and not current_has_sequence:
                    raise ValueError(
                        f"{source}:{line_number - 1}: FASTA record "
                        f"{current_id!r} has no sequence"
                    )
                header = text[1:]
                identifier = header.split(maxsplit=1)[0] if header else ""
                if not identifier:
                    raise ValueError(f"{source}:{line_number}: empty FASTA identifier")
                if identifier in identifiers:
                    raise ValueError(
                        f"{source}:{line_number}: duplicate FASTA identifier "
                        f"{identifier!r}"
                    )
                identifiers.add(identifier)
                current_id = identifier
                current_has_sequence = False
                count += 1
            elif current_id is None:
                raise ValueError(
                    f"{source}:{line_number}: sequence appears before a FASTA header"
                )
            elif text.strip():
                current_has_sequence = True
    if current_id is None:
        raise ValueError(f"{source}: FASTA contains no records")
    if not current_has_sequence:
        raise ValueError(f"{source}: FASTA record {current_id!r} has no sequence")
    return count


def fasta_lines(records: Iterable[FastaRecord]) -> Tuple[str, ...]:
    """Return deterministic two-line FASTA output preserving source headers."""

    lines: list[str] = []
    seen = set()
    for record in records:
        if record.identifier in seen:
            raise ValueError(f"duplicate FASTA output identifier {record.identifier!r}")
        seen.add(record.identifier)
        lines.extend((f">{record.header}", record.sequence))
    return tuple(lines)


def _make_record(
    header: str,
    sequence_parts: list[str],
    identifiers: set[str],
    source: str,
    line_number: int,
) -> FastaRecord:
    identifier = header.split(maxsplit=1)[0] if header else ""
    if not identifier:
        raise ValueError(f"{source}:{line_number}: empty FASTA identifier")
    if identifier in identifiers:
        raise ValueError(
            f"{source}:{line_number}: duplicate FASTA identifier {identifier!r}"
        )
    identifiers.add(identifier)
    sequence = "".join(sequence_parts)
    if not sequence:
        raise ValueError(
            f"{source}:{line_number}: FASTA record {identifier!r} has no sequence"
        )
    return FastaRecord(identifier, header, sequence)
