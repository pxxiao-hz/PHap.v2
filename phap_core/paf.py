"""Strict, lossless parsing of Pairwise mApping Format (PAF) records.

PAF coordinates are represented exactly as specified by the format: 0-based
half-open intervals on the forward strands of the query and target. Optional
fields are parsed by their two-character tag names and retain their original
order for lossless output.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Tuple, Union


class PafFormatError(ValueError):
    """Raised when a PAF record violates the format contract."""


@dataclass(frozen=True)
class PafTag:
    """One PAF optional field."""

    name: str
    type_code: str
    value: str
    raw: str


@dataclass(frozen=True)
class PafRecord:
    """One validated PAF record with its original fields."""

    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    matching_bases: int
    alignment_block_length: int
    mapping_quality: int
    tags: Tuple[PafTag, ...]
    fields: Tuple[str, ...]
    source: str
    line_number: int

    def tag(self, name: str) -> Optional[PafTag]:
        """Return an optional field by tag name."""

        return next((tag for tag in self.tags if tag.name == name), None)

    @property
    def alignment_type(self) -> Optional[str]:
        """Return the minimap2 ``tp:A`` value, or ``None`` when absent."""

        tag = self.tag("tp")
        return tag.value if tag is not None else None

    @property
    def is_primary(self) -> bool:
        """Whether the record explicitly declares ``tp:A:P``."""

        return self.alignment_type == "P"

    def to_line(self) -> str:
        """Return the original record without its line ending."""

        return "\t".join(self.fields)


@dataclass(frozen=True)
class PafRecordDecision:
    """Auditable primary-alignment filtering decision."""

    source: str
    line_number: int
    query_name: str
    target_name: str
    alignment_type: Optional[str]
    status: str
    reason: str


def parse_paf_line(
    line: str,
    *,
    source: str = "<memory>",
    line_number: int = 1,
) -> PafRecord:
    """Parse and validate one non-empty PAF line."""

    text = line.rstrip("\r\n")
    fields = tuple(text.split("\t"))
    if len(fields) < 12:
        raise _error(source, line_number, "expected at least 12 tab-delimited fields")
    if not fields[0] or not fields[5]:
        raise _error(source, line_number, "query and target names must not be empty")

    query_length = _integer(fields[1], source, line_number, "query length")
    query_start = _integer(fields[2], source, line_number, "query start")
    query_end = _integer(fields[3], source, line_number, "query end")
    strand = fields[4]
    target_length = _integer(fields[6], source, line_number, "target length")
    target_start = _integer(fields[7], source, line_number, "target start")
    target_end = _integer(fields[8], source, line_number, "target end")
    matching_bases = _integer(fields[9], source, line_number, "matching bases")
    alignment_block_length = _integer(
        fields[10],
        source,
        line_number,
        "alignment block length",
    )
    mapping_quality = _integer(fields[11], source, line_number, "mapping quality")

    if query_length < 1 or target_length < 1:
        raise _error(source, line_number, "query and target lengths must be positive")
    _validate_interval(
        query_start,
        query_end,
        query_length,
        source,
        line_number,
        "query",
    )
    _validate_interval(
        target_start,
        target_end,
        target_length,
        source,
        line_number,
        "target",
    )
    if strand not in {"+", "-"}:
        raise _error(source, line_number, "strand must be '+' or '-'")
    if alignment_block_length < 0:
        raise _error(source, line_number, "alignment block length must be non-negative")
    if not 0 <= matching_bases <= alignment_block_length:
        raise _error(
            source,
            line_number,
            "matching bases must be between zero and alignment block length",
        )
    if not 0 <= mapping_quality <= 255:
        raise _error(source, line_number, "mapping quality must be in [0, 255]")

    tags = _parse_tags(fields[12:], source, line_number)
    tp_tag = next((tag for tag in tags if tag.name == "tp"), None)
    if tp_tag is not None and (
        tp_tag.type_code != "A" or len(tp_tag.value) != 1
    ):
        raise _error(source, line_number, "tp tag must have form tp:A:<character>")

    return PafRecord(
        query_name=fields[0],
        query_length=query_length,
        query_start=query_start,
        query_end=query_end,
        strand=strand,
        target_name=fields[5],
        target_length=target_length,
        target_start=target_start,
        target_end=target_end,
        matching_bases=matching_bases,
        alignment_block_length=alignment_block_length,
        mapping_quality=mapping_quality,
        tags=tags,
        fields=fields,
        source=source,
        line_number=line_number,
    )


def parse_paf_lines(
    lines: Iterable[str],
    *,
    source: str = "<memory>",
) -> Tuple[PafRecord, ...]:
    """Parse all non-empty, non-comment PAF lines or raise without a result."""

    records = []
    for line_number, line in enumerate(lines, start=1):
        if not line.strip() or line.startswith("#"):
            continue
        records.append(
            parse_paf_line(
                line,
                source=source,
                line_number=line_number,
            )
        )
    return tuple(records)


def read_paf(path: Union[str, Path]) -> Tuple[PafRecord, ...]:
    """Read a PAF file with path and line numbers in format errors."""

    source = str(path)
    with Path(path).open(encoding="utf-8") as handle:
        return parse_paf_lines(handle, source=source)


def select_primary_records(
    records: Iterable[PafRecord],
    *,
    require_tp: bool,
) -> Tuple[Tuple[PafRecord, ...], Tuple[PafRecordDecision, ...]]:
    """Apply an explicit minimap2 primary-alignment policy.

    With ``require_tp=True``, records missing ``tp`` are excluded rather than
    guessed to be primary. Secondary and other non-primary ``tp`` values are
    always excluded.
    """

    accepted = []
    decisions = []
    for record in records:
        alignment_type = record.alignment_type
        if alignment_type == "P":
            status = "accepted"
            reason = "primary_alignment"
            accepted.append(record)
        elif alignment_type is None and not require_tp:
            status = "accepted"
            reason = "tp_not_required"
            accepted.append(record)
        elif alignment_type is None:
            status = "excluded"
            reason = "missing_tp_tag"
        elif alignment_type == "S":
            status = "excluded"
            reason = "secondary_alignment"
        else:
            status = "excluded"
            reason = "non_primary_alignment"
        decisions.append(
            PafRecordDecision(
                source=record.source,
                line_number=record.line_number,
                query_name=record.query_name,
                target_name=record.target_name,
                alignment_type=alignment_type,
                status=status,
                reason=reason,
            )
        )
    return (
        tuple(sorted(accepted, key=_record_sort_key)),
        tuple(sorted(decisions, key=lambda row: (row.source, row.line_number))),
    )


def _record_sort_key(record: PafRecord) -> Tuple[object, ...]:
    return (
        record.query_name,
        record.target_name,
        record.query_start,
        record.query_end,
        record.strand,
        record.target_start,
        record.target_end,
        record.line_number,
    )


def _integer(
    value: str,
    source: str,
    line_number: int,
    field_name: str,
) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise _error(source, line_number, f"{field_name} must be an integer") from exc


def _validate_interval(
    start: int,
    end: int,
    length: int,
    source: str,
    line_number: int,
    label: str,
) -> None:
    if start < 0 or end < start or end > length:
        raise _error(
            source,
            line_number,
            f"invalid 0-based half-open {label} interval [{start}, {end}) "
            f"for length {length}",
        )


def _parse_tags(
    fields: Tuple[str, ...],
    source: str,
    line_number: int,
) -> Tuple[PafTag, ...]:
    tags = []
    seen = set()
    for field in fields:
        parts = field.split(":", 2)
        if (
            len(parts) != 3
            or len(parts[0]) != 2
            or len(parts[1]) != 1
            or not parts[2]
        ):
            raise _error(source, line_number, f"malformed optional field {field!r}")
        name, type_code, value = parts
        if name in seen:
            raise _error(source, line_number, f"duplicate optional tag {name!r}")
        seen.add(name)
        tags.append(PafTag(name, type_code, value, field))
    return tuple(tags)


def _error(source: str, line_number: int, reason: str) -> PafFormatError:
    return PafFormatError(f"{source}:{line_number}: {reason}")
