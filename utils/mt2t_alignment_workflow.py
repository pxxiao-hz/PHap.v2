"""Build pairwise mT2T PAF evidence without stale caches or shell interpolation."""

from __future__ import annotations

import hashlib
import shlex
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Iterator, Sequence, Tuple

from phap_core import __version__
from phap_core.atomic_io import atomic_write_lines
from phap_core.fasta import FastaRecord, read_fasta
from phap_core.runner import run_command


def build_pairwise_overlap_paf(
    *,
    fasta_path: str,
    output_directory: str,
    min_contig_length: int,
    min_distance: float,
    threads_per_alignment: int,
    max_processes: int,
    cpu_budget: int,
    tool_versions: Sequence[Tuple[str, str]],
) -> Path:
    """Always rebuild Mash/minimap2 evidence in an isolated staging directory."""

    if min_contig_length < 1:
        raise ValueError("min_contig_length must be positive")
    if not 0.0 <= min_distance <= 1.0:
        raise ValueError("min_distance must be in [0, 1]")
    if threads_per_alignment < 1 or max_processes < 1 or cpu_budget < 1:
        raise ValueError("alignment threads, process count, and CPU budget must be positive")
    if threads_per_alignment * max_processes > cpu_budget:
        raise ValueError("alignment concurrency exceeds the CPU budget")
    records = tuple(
        record
        for record in read_fasta(fasta_path)
        if len(record.sequence) >= min_contig_length
    )
    output_root = Path(output_directory)
    alignment_root = output_root / "03.alignment"
    alignment_root.mkdir(parents=True, exist_ok=True)
    merge_paf = alignment_root / "merge.paf"
    upstream_manifest = output_root / "mt2t_upstream_manifest.tsv"
    upstream_manifest.unlink(missing_ok=True)
    commands: list[str] = []
    similar_pairs: list[Tuple[int, int, float]] = []

    if len(records) < 2:
        atomic_write_lines(merge_paf, ())
    else:
        with tempfile.TemporaryDirectory(
            prefix=".mt2t-alignment-",
            dir=output_root,
        ) as temporary_directory:
            staging = Path(temporary_directory)
            fasta_files = _write_safe_split_fastas(records, staging)
            sketch_files = tuple(
                staging / f"sketch_{index:06d}.msh"
                for index in range(1, len(records) + 1)
            )
            sketch_commands = [
                (
                    "mash",
                    "sketch",
                    "-o",
                    str(sketch.with_suffix("")),
                    str(fasta_file),
                )
                for sketch, fasta_file in zip(sketch_files, fasta_files)
            ]
            commands.extend(
                _portable_command(command, staging) for command in sketch_commands
            )
            _run_parallel(sketch_commands, max_processes=max_processes)

            pair_rows = [
                (left, right)
                for left in range(len(records))
                for right in range(left + 1, len(records))
            ]

            def mash_distance(pair: Tuple[int, int]) -> Tuple[int, int, float]:
                left, right = pair
                command = (
                    "mash",
                    "dist",
                    str(sketch_files[left]),
                    str(sketch_files[right]),
                )
                completed = run_command(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                output = completed.stdout
                text = (
                    output.decode("utf-8", errors="strict")
                    if isinstance(output, bytes)
                    else str(output)
                )
                fields = text.strip().split("\t")
                if len(fields) < 3:
                    raise ValueError("mash dist returned a malformed row")
                return left, right, float(fields[2])

            commands.extend(
                _portable_command(
                    (
                        "mash",
                        "dist",
                        str(sketch_files[left]),
                        str(sketch_files[right]),
                    ),
                    staging,
                )
                for left, right in pair_rows
            )
            with ThreadPoolExecutor(
                max_workers=min(max_processes, len(pair_rows))
            ) as executor:
                distance_rows = tuple(executor.map(mash_distance, pair_rows))
            similar_pairs = sorted(
                row for row in distance_rows if row[2] < min_distance
            )

            pair_pafs: list[Path] = []

            def align_pair(row: Tuple[int, int, float]) -> Path:
                left, right, _ = row
                output_path = staging / f"pair_{left:06d}_{right:06d}.paf"
                command = (
                    "minimap2",
                    "-cx",
                    "asm5",
                    "-t",
                    str(threads_per_alignment),
                    str(fasta_files[left]),
                    str(fasta_files[right]),
                )
                with output_path.open("w", encoding="utf-8", newline="\n") as output:
                    run_command(command, stdout=output)
                return output_path

            commands.extend(
                _portable_command(
                    (
                        "minimap2",
                        "-cx",
                        "asm5",
                        "-t",
                        str(threads_per_alignment),
                        str(fasta_files[left]),
                        str(fasta_files[right]),
                    ),
                    staging,
                )
                for left, right, _ in similar_pairs
            )
            if similar_pairs:
                with ThreadPoolExecutor(
                    max_workers=min(max_processes, len(similar_pairs))
                ) as executor:
                    pair_pafs = list(executor.map(align_pair, similar_pairs))
            atomic_write_lines(merge_paf, _paf_lines(pair_pafs))

    atomic_write_lines(
        output_root / "mt2t_split_id_map.tsv",
        (
            "safe_file_ID\tsource_ID\tsource_header",
            *(
                f"contig_{index:06d}\t{record.identifier}\t{record.header}"
                for index, record in enumerate(records, start=1)
            ),
        ),
    )
    atomic_write_lines(
        upstream_manifest,
        (
            "key\tvalue",
            f"phap_version\t{__version__}",
            f"fasta\t{fasta_path}",
            f"fasta_sha256\t{_sha256(Path(fasta_path))}",
            f"min_contig_length\t{min_contig_length}",
            f"min_distance\t{min_distance:.12g}",
            f"threads_per_alignment\t{threads_per_alignment}",
            f"max_processes\t{max_processes}",
            f"cpu_budget\t{cpu_budget}",
            f"eligible_contigs\t{len(records)}",
            f"similar_pairs\t{len(similar_pairs)}",
            *(f"tool_{name}\t{version}" for name, version in sorted(tool_versions)),
            *(f"command_{index:06d}\t{command}" for index, command in enumerate(commands, 1)),
            f"merge_paf_sha256\t{_sha256(merge_paf)}",
        ),
    )
    return merge_paf


def _write_safe_split_fastas(
    records: Sequence[FastaRecord],
    directory: Path,
) -> Tuple[Path, ...]:
    paths = []
    for index, record in enumerate(records, start=1):
        path = directory / f"contig_{index:06d}.fa"
        atomic_write_lines(path, (f">{record.header}", record.sequence))
        paths.append(path)
    return tuple(paths)


def _run_parallel(
    commands: Sequence[Sequence[str]],
    *,
    max_processes: int,
) -> None:
    if not commands:
        return
    with ThreadPoolExecutor(max_workers=min(max_processes, len(commands))) as executor:
        tuple(executor.map(run_command, commands))


def _paf_lines(paths: Sequence[Path]) -> Iterator[str]:
    for path in sorted(paths):
        with path.open(encoding="utf-8") as source:
            for line in source:
                if line.strip():
                    yield line.rstrip("\r\n")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _portable_command(command: Sequence[str], staging: Path) -> str:
    return shlex.join(command).replace(str(staging), "$STAGING")
