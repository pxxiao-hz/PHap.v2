"""Shared validation and external-command execution helpers."""

from __future__ import annotations

import logging
import os
import shlex
import shutil
import subprocess
from collections.abc import Iterable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import IO, Any, Optional, Union

PathLike = Union[os.PathLike[str], str]
Stream = Optional[Union[IO[Any], int]]

LOGGER = logging.getLogger("phap.commands")


class PreflightError(RuntimeError):
    """Raised when required inputs or executables are unavailable."""


def require_input_files(paths: Iterable[PathLike]) -> None:
    """Require all input paths to be existing regular files."""
    missing = [str(path) for path in paths if not Path(path).is_file()]
    if missing:
        raise PreflightError("missing input file(s): " + ", ".join(missing))


def require_tools(tools: Iterable[str]) -> None:
    """Require external executables to be discoverable through PATH."""
    missing = sorted(tool for tool in set(tools) if shutil.which(tool) is None)
    if missing:
        raise PreflightError(
            "required executable(s) not found on PATH: " + ", ".join(missing)
        )


def run_command(
    args: Sequence[PathLike],
    *,
    cwd: Optional[PathLike] = None,
    env: Optional[Mapping[str, str]] = None,
    stdout: Stream = None,
    stderr: Stream = None,
    check: bool = True,
) -> subprocess.CompletedProcess[Any]:
    """Run one command without a shell and propagate failures by default."""
    command = [str(arg) for arg in args]
    LOGGER.info("Running: %s", shlex.join(command))
    return subprocess.run(
        command,
        cwd=cwd,
        env=None if env is None else dict(env),
        stdout=stdout,
        stderr=stderr,
        check=check,
    )


def run_shell_command(
    command: str,
    *,
    cwd: Optional[PathLike] = None,
    env: Optional[Mapping[str, str]] = None,
) -> None:
    """Run a legacy shell pipeline while preserving its failure status."""
    LOGGER.info("Running shell pipeline: %s", command)
    subprocess.run(
        command,
        shell=True,
        executable="/bin/bash",
        cwd=cwd,
        env=None if env is None else dict(env),
        check=True,
    )


def run_shell_commands_parallel(
    commands: Iterable[str],
    max_workers: int,
    env: Optional[Mapping[str, str]] = None,
) -> None:
    """Run shell pipelines concurrently and re-raise every worker failure."""
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")

    command_list = list(commands)
    if not command_list:
        return

    with ThreadPoolExecutor(max_workers=min(max_workers, len(command_list))) as executor:
        futures = [
            executor.submit(run_shell_command, command, env=env)
            for command in command_list
        ]
        for future in futures:
            future.result()
