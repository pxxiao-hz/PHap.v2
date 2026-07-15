"""Content-addressed minimap2 alignment stage for PHap clustering."""

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, Optional

from phap_core import __version__
from phap_core.atomic_io import atomic_text_writer
from phap_core.runner import run_command
from phap_core.stage_manifest import (
    CacheValidation,
    StageInput,
    StageSpec,
    capture_stage,
    invalidate_stage_manifest,
    stage_lock,
    validate_stage_cache,
    write_stage_manifest,
)


CommandRunner = Callable[..., Any]
ALIGNMENT_IMPLEMENTATION = "cluster-locus-alignment-v1"


def build_or_reuse_locus_alignment(
    *,
    mt2t_fasta: str,
    unitig_fasta: str,
    output_paf: str,
    manifest_path: str,
    threads: int,
    minimap2_version: str,
    minimap2_executable: Optional[str] = None,
    phap_version: str = __version__,
    runner: CommandRunner = run_command,
) -> CacheValidation:
    """Build an atomic PAF or reuse it only after full manifest validation."""

    if threads < 1:
        raise ValueError("threads must be positive")
    _require_distinct_paths(
        {
            "mt2t_fasta": mt2t_fasta,
            "unitig_fasta": unitig_fasta,
            "output_paf": output_paf,
            "manifest": manifest_path,
        }
    )
    executable = _resolve_executable(minimap2_executable)
    if not minimap2_version.strip():
        raise ValueError("minimap2_version must not be empty")
    command = (
        executable,
        "-cx",
        "asm5",
        "-t",
        str(threads),
        str(Path(mt2t_fasta).resolve()),
        str(Path(unitig_fasta).resolve()),
    )
    spec = StageSpec(
        stage="cluster.putg_vs_mt2t_alignment",
        inputs=(
            StageInput("mt2t_fasta", mt2t_fasta),
            StageInput("unitig_fasta", unitig_fasta),
        ),
        parameters=(
            ("command", _command_text(command)),
            ("implementation", ALIGNMENT_IMPLEMENTATION),
            ("preset", "asm5"),
            ("threads", str(threads)),
            ("tool_executable", executable),
        ),
        tool_versions=(("minimap2", minimap2_version.strip()),),
        phap_version=phap_version,
    )
    outputs = {"paf": output_paf}

    with stage_lock(manifest_path):
        snapshot = capture_stage(spec)
        validation = validate_stage_cache(snapshot, manifest_path, outputs)
        if validation.hit:
            return validation
        invalidate_stage_manifest(manifest_path)
        with atomic_text_writer(output_paf) as output:
            runner(command, stdout=output)
        write_stage_manifest(snapshot, manifest_path, outputs)
        return validation


def minimap2_version(
    executable: Optional[str] = None,
    *,
    runner: CommandRunner = run_command,
) -> tuple[str, str]:
    """Return the resolved executable and its first version-output line."""

    resolved = _resolve_executable(executable)
    completed = runner(
        [resolved, "--version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    output = completed.stdout
    text = (
        output.decode("utf-8", errors="replace")
        if isinstance(output, bytes)
        else str(output)
    )
    lines = text.strip().splitlines()
    if not lines or not lines[0]:
        raise ValueError("minimap2 --version returned no version text")
    return resolved, lines[0]


def _resolve_executable(executable: Optional[str]) -> str:
    candidate = executable or "minimap2"
    resolved = shutil.which(candidate)
    if resolved is None:
        raise ValueError(f"minimap2 executable not found on PATH: {candidate}")
    return str(Path(resolved).resolve())


def _command_text(command: Sequence[str]) -> str:
    return "\x1f".join(command)


def _require_distinct_paths(paths: dict[str, str]) -> None:
    by_path: dict[Path, str] = {}
    for role, raw_path in paths.items():
        path = Path(raw_path).resolve()
        previous_role = by_path.get(path)
        if previous_role is not None:
            raise ValueError(
                f"stage paths must be distinct: {previous_role} and {role} "
                f"both resolve to {path}"
            )
        by_path[path] = role
