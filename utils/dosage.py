#!/usr/bin/env python3
"""Shared parsing and naming helpers for copy-number dosage classes."""

from __future__ import annotations

import re
from typing import Optional


# Keep the historical labels in outputs where a conventional name exists.
_LABEL_BY_DOSAGE = {
    1: "haplotig",
    2: "diplotig",
    3: "triplotig",
    4: "tetraplotig",
    5: "pentaplotig",
    6: "hexaplotig",
}
DOSAGE_BY_TYPE = {label: dosage for dosage, label in _LABEL_BY_DOSAGE.items()}

_GENERIC_PATTERNS = (
    re.compile(r"^([1-9][0-9]*)$"),
    re.compile(r"^dosage[_-]?([1-9][0-9]*)$"),
    re.compile(r"^copy[_-]?([1-9][0-9]*)$"),
    re.compile(r"^([1-9][0-9]*)x$"),
)


def dosage_from_contig_type(contig_type: object) -> Optional[int]:
    """Return a positive integer dosage for a supported contig-type label.

    Historical labels remain accepted.  A bare integer, ``dosage_N``,
    ``copy_N``, and ``Nx``
    make the file format usable at arbitrary ploidy without adding another
    hard-coded biological name for every copy number.
    """

    if contig_type is None:
        return None
    label = str(contig_type).strip().lower()
    if label in DOSAGE_BY_TYPE:
        return DOSAGE_BY_TYPE[label]
    for pattern in _GENERIC_PATTERNS:
        match = pattern.fullmatch(label)
        if match:
            return int(match.group(1))
    return None


def contig_type_for_dosage(dosage: int) -> str:
    """Return a backward-compatible label for a positive integer dosage."""

    if dosage < 1:
        raise ValueError("dosage must be positive")
    return _LABEL_BY_DOSAGE.get(dosage, f"dosage_{dosage}")


def validate_dosage(contig_type: object, ploidy: int) -> Optional[int]:
    """Parse a dosage label and return it only when it fits the ploidy."""

    dosage = dosage_from_contig_type(contig_type)
    if dosage is None or dosage > ploidy:
        return None
    return dosage
