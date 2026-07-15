"""Atomic text-output helpers shared by auditable workflow stages."""

from __future__ import annotations

import os
from collections.abc import Iterable
from pathlib import Path
from typing import Union

PathLike = Union[str, os.PathLike[str]]


def atomic_write_lines(path: PathLike, lines: Iterable[str]) -> None:
    """Write newline-terminated text and atomically replace the destination."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as output:
            for line in lines:
                output.write(line)
                output.write("\n")
        os.replace(temporary, target)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
