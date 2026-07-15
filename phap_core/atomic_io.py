"""Atomic text-output helpers shared by auditable workflow stages."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Iterable
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO, Union

PathLike = Union[str, os.PathLike[str]]


def atomic_write_lines(path: PathLike, lines: Iterable[str]) -> None:
    """Write newline-terminated text and atomically replace the destination."""

    with atomic_text_writer(path) as output:
        for line in lines:
            output.write(line)
            output.write("\n")


@contextmanager
def atomic_text_writer(path: PathLike) -> Iterator[TextIO]:
    """Yield a temporary text handle and replace ``path`` only on success."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=target.parent,
        prefix=f".{target.name}.",
        suffix=".tmp",
        text=True,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as output:
            yield output
        os.replace(temporary, target)
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        temporary.unlink(missing_ok=True)
        raise
