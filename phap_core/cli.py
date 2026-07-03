"""Top-level PHap command dispatcher."""

from __future__ import annotations

import subprocess
import sys
from collections.abc import Sequence
from typing import Optional

from . import __update_time__, __version__


COMMANDS = {
    "dosage": (
        "scripts.phap_dosage",
        "Infer unitig copy dosage from window-level read depth.",
    ),
    "mt2t": (
        "scripts.phap_mT2T",
        "Generate a mosaic T2T reference from a primary contig assembly.",
    ),
    "cluster": (
        "scripts.phap_cluster",
        "Build the allelic table and cluster unitigs with Hi-C evidence.",
    ),
    "phase_reads": (
        "scripts.phap_phase_reads",
        "Phase reads, assemble each haplotype, and scaffold the assemblies.",
    ),
}


def format_help() -> str:
    """Return stable top-level CLI help without importing scientific modules."""
    command_lines = "\n".join(
        f"  {name:<12}{description}" for name, (_, description) in COMMANDS.items()
    )
    return (
        "Usage: phap <command> [options]\n\n"
        "Commands:\n"
        f"{command_lines}\n"
        "  version     Print version information.\n\n"
        "Run 'phap <command> --help' for command-specific options."
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Dispatch a PHap subcommand and return its exit status."""
    args = list(sys.argv[1:] if argv is None else argv)

    if not args or args[0] in {"-h", "--help"}:
        print(format_help())
        return 0

    if args[0] in {"version", "-v", "--version", "-version"}:
        print(f"Version: {__version__}\nUpdate Time: {__update_time__}")
        return 0

    command = COMMANDS.get(args[0])
    if command is None:
        print(f"phap: error: unknown command '{args[0]}'", file=sys.stderr)
        print(format_help(), file=sys.stderr)
        return 2

    module_name = command[0]
    completed = subprocess.run(
        [sys.executable, "-m", module_name, *args[1:]],
        check=False,
    )
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
