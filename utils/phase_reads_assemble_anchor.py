#!/usr/bin/env python3
"""Phase reads, assemble groups, and scaffold them with auditable stages."""

from __future__ import annotations

import argparse
import errno
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from phase_reads_assignment import (
    AssignmentParameters,
    assign_hic_fast_v1,
    assign_hic_pairs_disk,
    assign_long_reads_disk,
    collect_evidence_disk,
    dispatch_paired_fastq_fast_v1,
    dispatch_paired_fastq_database,
    dispatch_single_fastq_database,
    dispatch_single_fastq_seqkit,
    group_sort_key,
    group_target_lengths,
    parse_contig_types,
    parse_group_file,
    validate_group_dosage,
    write_assignments_database,
    write_json,
)


STAGES = ("assign", "extract", "assemble", "scaffold")
LOGGER = logging.getLogger("phap.phase_reads")


class ProgressLog:
    def __init__(self, label, every):
        self.label = label
        self.every = max(every, 1)
        self.last = 0
        self.lock = threading.Lock()

    def __call__(self, processed, total, started, detail=None):
        with self.lock:
            if processed < self.last + self.every and not (
                total is not None and processed >= total
            ):
                return
            self.last = processed
            elapsed = max(time.monotonic() - started, 0.001)
            rate = processed / elapsed
            suffix = f", {detail}" if detail else ""
            if total:
                percent = 100.0 * processed / total
                LOGGER.info(
                    "%s: %s/%s (%.1f%%), %.1f records/s%s",
                    self.label, f"{processed:,}", f"{total:,}", percent, rate,
                    suffix,
                )
            else:
                LOGGER.info(
                    "%s: %s processed, %.1f records/s%s",
                    self.label, f"{processed:,}", rate, suffix,
                )


def add_path_argument(parser, hyphen_name, underscore_name=None, **kwargs):
    names = [hyphen_name]
    if underscore_name:
        names.append(underscore_name)
    parser.add_argument(*names, **kwargs)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Auditable haplotype read assignment, assembly, and scaffolding"
    )
    add_path_argument(parser, "--bam-hifi", "--bam_hifi", type=Path)
    add_path_argument(parser, "--bam-hic", "--bam_hic", type=Path)
    add_path_argument(parser, "--bam-ont", "--bam_ont", type=Path)
    add_path_argument(
        parser, "--contig-type", "--contig_type", type=Path
    )
    parser.add_argument(
        "--group",
        type=Path,
        default=Path("02.cluster.v2/05.rescue/group.reassignment.cluster.txt"),
    )
    parser.add_argument("--hifi", type=Path)
    parser.add_argument("--ont", type=Path)
    parser.add_argument("--hic1", type=Path)
    parser.add_argument("--hic2", type=Path)
    add_path_argument(
        parser,
        "--output-dir",
        "--output_dir",
        type=Path,
        default=Path("03.phase_reads"),
    )
    add_path_argument(
        parser,
        "--temp-dir",
        "--temp_dir",
        type=Path,
        help=(
            "Directory for PHap working files and SQLite/external-tool temporary "
            "files; final outputs remain under --output-dir"
        ),
    )
    parser.add_argument(
        "--scaffold-only",
        action="store_true",
        help=(
            "Run only HapHiC scaffolding from existing per-group assemblies "
            "and extracted Hi-C FASTQs; requires --assembly-dir and "
            "--hic-reads-dir"
        ),
    )
    add_path_argument(
        parser,
        "--assembly-dir",
        "--assembly_dir",
        type=Path,
        help=(
            "Existing per-group assembly directory containing "
            "GROUP/GROUP.asm.bp.p_ctg.gfa (scaffold-only mode)"
        ),
    )
    add_path_argument(
        parser,
        "--hic-reads-dir",
        "--hic_reads_dir",
        type=Path,
        help=(
            "Existing directory containing GROUP.Hi-C.1.fq.gz and "
            "GROUP.Hi-C.2.fq.gz (scaffold-only mode)"
        ),
    )
    parser.add_argument("--stop-after", choices=STAGES, default="scaffold")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--rerun-from", choices=STAGES,
        help="With --resume, keep earlier checkpoints and rerun from this stage",
    )
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--groups", nargs="+",
        help="Output only these groups; sibling groups are still used for fair assignment",
    )
    selection.add_argument(
        "--chromosomes", nargs="+",
        help="Run all haplotype groups for these chromosomes",
    )
    parser.add_argument(
        "--data-types", "--data_types", nargs="+",
        choices=("hifi", "ont", "hic"),
        default=["hifi", "ont", "hic"],
        help="Sequencing types to process [hifi ont hic]",
    )
    parser.add_argument("--seed", type=int, default=100)
    parser.add_argument(
        "--progress-every", type=int, default=1_000_000,
        help="Report record progress at this interval [1000000]",
    )
    parser.add_argument(
        "--collapsed-policy", choices=("balanced", "defer"), default="balanced"
    )
    parser.add_argument(
        "--unknown-contig-type-policy", choices=("group", "error"), default="group"
    )
    parser.add_argument("--min-group-margin", type=float, default=0.10)
    parser.add_argument("--balanced-score-tolerance", type=float, default=0.02)
    parser.add_argument("--hifi-min-alignment-length", type=int, default=1000)
    parser.add_argument("--hifi-min-identity", type=float, default=0.95)
    parser.add_argument("--ont-min-alignment-length", type=int, default=1000)
    parser.add_argument("--ont-min-identity", type=float, default=0.75)
    parser.add_argument("--hic-min-alignment-length", type=int, default=50)
    parser.add_argument("--hic-min-identity", type=float, default=0.90)
    parser.add_argument(
        "--hic-assignment-backend",
        choices=("fast-v1", "sqlite"),
        default="fast-v1",
        help=(
            "Hi-C assignment engine: corrected v1 primary-read1/RNEXT streaming "
            "without evidence SQLite, or the detailed pair-evidence SQLite "
            "engine [fast-v1]"
        ),
    )
    add_path_argument(parser, "--ont-length", "--ont_length", type=int, default=1)
    add_path_argument(parser, "--ont-quality", "--ont_quality", type=float, default=0.0)
    parser.add_argument(
        "--extract-backend", choices=("auto", "python", "seqkit"), default="auto",
        help=(
            "FASTQ extraction backend for HiFi/ONT: auto uses seqkit/gawk/pigz "
            "when available and otherwise Python [auto]"
        ),
    )
    parser.add_argument(
        "--extract-threads", type=int, default=8,
        help="Threads used by seqkit input decompression [8]",
    )
    parser.add_argument(
        "--jobs", "--process", dest="jobs", type=int, default=4,
        help="Maximum concurrent group jobs [4]",
    )
    parser.add_argument(
        "--threads-per-job", "--threads", dest="threads_per_job", type=int, default=10,
        help="Threads supplied to each external group job [10]",
    )
    parser.add_argument("--haphic-processes", type=int, default=1)
    parser.add_argument("--haphic-nx", type=int, default=100)
    parser.add_argument("--haphic-nm", type=int, default=3)
    parser.add_argument("--run-juicebox", action="store_true")
    parser.add_argument("--hifiasm", default="hifiasm")
    parser.add_argument("--bwa", default="bwa")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--samblaster", default="samblaster")
    parser.add_argument("--filter-bam", default="filter_bam")
    parser.add_argument("--haphic", default="haphic")
    parser.add_argument("--seqkit", default="seqkit")
    parser.add_argument("--pigz", default="pigz")
    parser.add_argument("--gawk", default="gawk")
    return parser.parse_args(argv)


def resolved_inputs(args):
    for name in (
        "bam_hifi", "bam_ont", "bam_hic", "contig_type", "group",
        "hifi", "ont", "hic1", "hic2", "output_dir", "temp_dir",
        "assembly_dir", "hic_reads_dir",
    ):
        value = getattr(args, name)
        if value is not None:
            setattr(args, name, value.resolve())


def validate_args(args):
    if (
        args.jobs < 1 or args.threads_per_job < 1
        or args.haphic_processes < 1 or args.extract_threads < 1
    ):
        raise ValueError(
            "jobs, threads-per-job, haphic-processes, and extract-threads "
            "must be positive"
        )
    if args.progress_every < 1:
        raise ValueError("--progress-every must be positive")
    if args.rerun_from and not args.resume:
        raise ValueError("--rerun-from requires --resume")
    if args.rerun_from and STAGES.index(args.rerun_from) > STAGES.index(args.stop_after):
        raise ValueError("--rerun-from cannot be later than --stop-after")
    if not 0 <= args.min_group_margin <= 1:
        raise ValueError("--min-group-margin must be between zero and one")
    if not 0 <= args.balanced_score_tolerance <= 1:
        raise ValueError("--balanced-score-tolerance must be between zero and one")
    if args.ont_length < 0 or args.ont_quality < 0:
        raise ValueError("ONT length and quality filters cannot be negative")
    if args.scaffold_only:
        if args.stop_after != "scaffold":
            raise ValueError("--scaffold-only requires --stop-after scaffold")
        if args.rerun_from not in (None, "scaffold"):
            raise ValueError(
                "--scaffold-only supports only --rerun-from scaffold"
            )
        missing_arguments = [
            label for value, label in (
                (args.assembly_dir, "--assembly-dir"),
                (args.hic_reads_dir, "--hic-reads-dir"),
            ) if value is None
        ]
        if missing_arguments:
            raise ValueError(
                "--scaffold-only requires " + " and ".join(missing_arguments)
            )
        required = [args.group, args.assembly_dir, args.hic_reads_dir]
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Missing scaffold-only inputs: {', '.join(missing)}")
        non_directories = [
            str(path) for path in (args.assembly_dir, args.hic_reads_dir)
            if not path.is_dir()
        ]
        if non_directories:
            raise ValueError(
                "Scaffold-only input is not a directory: "
                + ", ".join(non_directories)
            )
        if args.temp_dir is not None:
            if args.temp_dir.exists() and not args.temp_dir.is_dir():
                raise ValueError(f"--temp-dir is not a directory: {args.temp_dir}")
            args.temp_dir.mkdir(parents=True, exist_ok=True)
        return
    if args.assembly_dir is not None or args.hic_reads_dir is not None:
        raise ValueError(
            "--assembly-dir and --hic-reads-dir require --scaffold-only"
        )
    args.data_types = list(dict.fromkeys(args.data_types))
    if args.contig_type is None:
        raise ValueError("--contig-type is required unless --scaffold-only is used")
    required = [args.contig_type, args.group]
    for kind in args.data_types:
        value = getattr(args, f"bam_{kind}")
        if value is None:
            raise ValueError(f"--bam-{kind} is required when {kind} is selected")
        required.append(value)
    if STAGES.index(args.stop_after) >= STAGES.index("assemble") and "hifi" not in args.data_types:
        raise ValueError("hifi must be selected when running through assemble")
    if args.stop_after == "scaffold" and "hic" not in args.data_types:
        raise ValueError("hic must be selected when running through scaffold")
    if STAGES.index(args.stop_after) >= STAGES.index("extract"):
        raw_inputs = {
            "hifi": ((args.hifi, "--hifi"),),
            "ont": ((args.ont, "--ont"),),
            "hic": ((args.hic1, "--hic1"), (args.hic2, "--hic2")),
        }
        for value, label in (
            item for kind in args.data_types for item in raw_inputs[kind]
        ):
            if value is None:
                raise ValueError(f"{label} is required when running through extract")
            required.append(value)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing input files: {', '.join(missing)}")
    if args.temp_dir is not None:
        if args.temp_dir.exists() and not args.temp_dir.is_dir():
            raise ValueError(f"--temp-dir is not a directory: {args.temp_dir}")
        args.temp_dir.mkdir(parents=True, exist_ok=True)


def configure_temp_directory(args):
    """Confine Python, SQLite, and child-process temporary files when requested."""
    if args.temp_dir is None:
        return None
    path = args.temp_dir.resolve()
    os.environ["TMPDIR"] = str(path)
    os.environ["SQLITE_TMPDIR"] = str(path)
    # tempfile caches its selected directory after first use, so set it explicitly.
    tempfile.tempdir = str(path)
    token = hashlib.sha256(str(args.output_dir.resolve()).encode()).hexdigest()[:12]
    workspace = path / f"phap_phase_reads.{args.output_dir.name}.{token}"
    workspace.mkdir(parents=True, exist_ok=True)
    args._temp_workspace = workspace
    return workspace


def stage_temp_directory(args, stage):
    workspace = getattr(args, "_temp_workspace", None)
    if workspace is None:
        return None
    target = workspace / stage
    target.mkdir(parents=True, exist_ok=True)
    return target


def prune_empty_temp_workspace(args):
    workspace = getattr(args, "_temp_workspace", None)
    if workspace is None or not workspace.exists():
        return
    for path in sorted(
        (value for value in workspace.rglob("*") if value.is_dir()),
        key=lambda value: len(value.parts),
        reverse=True,
    ):
        try:
            path.rmdir()
        except OSError:
            pass
    try:
        workspace.rmdir()
    except OSError:
        pass


def _selection_values(values):
    if not values:
        return []
    result = []
    for value in values:
        result.extend(item for item in value.split(",") if item)
    return list(dict.fromkeys(result))


def subset_group_model(model, selected_groups):
    selected = set(selected_groups)
    groups = {
        group: members for group, members in model.groups.items() if group in selected
    }
    unitig_groups = {}
    for unitig, memberships in model.unitig_groups.items():
        retained = tuple(group for group in memberships if group in selected)
        if retained:
            unitig_groups[unitig] = retained
    return type(model)(
        groups=groups,
        unitig_groups=unitig_groups,
        group_chromosome={group: model.group_chromosome[group] for group in groups},
        group_number={group: model.group_number[group] for group in groups},
        ploidy=model.ploidy,
    )


def select_group_models(model, args):
    requested_groups = _selection_values(args.groups)
    requested_chromosomes = _selection_values(args.chromosomes)
    known_chromosomes = set(model.group_chromosome.values())
    if requested_groups:
        unknown = sorted(set(requested_groups) - set(model.groups))
        if unknown:
            raise ValueError(f"Unknown --groups values: {', '.join(unknown)}")
        output_groups = requested_groups
        chromosomes = {
            model.group_chromosome[group] for group in requested_groups
        }
    elif requested_chromosomes:
        unknown = sorted(set(requested_chromosomes) - known_chromosomes)
        if unknown:
            raise ValueError(
                f"Unknown --chromosomes values: {', '.join(unknown)}"
            )
        chromosomes = set(requested_chromosomes)
        output_groups = [
            group for group in model.groups
            if model.group_chromosome[group] in chromosomes
        ]
    else:
        chromosomes = known_chromosomes
        output_groups = list(model.groups)
    analysis_groups = [
        group for group in model.groups
        if model.group_chromosome[group] in chromosomes
    ]
    return (
        subset_group_model(model, analysis_groups),
        subset_group_model(model, output_groups),
    )


def resume_stage(args, stage):
    if not getattr(args, "resume", False):
        return False
    rerun_from = getattr(args, "rerun_from", None)
    if not rerun_from:
        return True
    return STAGES.index(stage) < STAGES.index(rerun_from)


def file_signature(path: Path):
    stat = path.stat()
    result = {
        "path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns
    }
    if stat.st_size <= 10 * 1024 * 1024:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        result["sha256"] = digest.hexdigest()
    return result


def scaffold_stage_configuration(args, inputs=None):
    return {
        "inputs": inputs or {},
        "parameters": {
            "bwa": args.bwa,
            "samtools": args.samtools,
            "samblaster": args.samblaster,
            "filter_bam": args.filter_bam,
            "haphic": args.haphic,
            "haphic_nx": args.haphic_nx,
            "haphic_nm": args.haphic_nm,
            "run_juicebox": args.run_juicebox,
        },
        "runtime": {
            "jobs": args.jobs,
            "threads_per_job": args.threads_per_job,
            "haphic_processes": args.haphic_processes,
        },
    }


def scaffold_only_group_inputs(args, groups):
    result = {}
    missing = []
    empty = []
    for group in sorted(groups, key=group_sort_key):
        paths = {
            "assembly_gfa": (
                args.assembly_dir / group / f"{group}.asm.bp.p_ctg.gfa"
            ),
            "hic1": args.hic_reads_dir / f"{group}.Hi-C.1.fq.gz",
            "hic2": args.hic_reads_dir / f"{group}.Hi-C.2.fq.gz",
        }
        for path in paths.values():
            if not path.exists():
                missing.append(str(path))
            elif not path.is_file() or path.stat().st_size == 0:
                empty.append(str(path))
        result[group] = paths
    if missing:
        raise FileNotFoundError(
            "Missing scaffold-only group inputs: " + ", ".join(missing)
        )
    if empty:
        raise ValueError(
            "Empty scaffold-only group inputs: " + ", ".join(empty)
        )
    return result


def scaffold_only_configuration(args, group_inputs):
    inputs = {
        "assembly_dir": str(args.assembly_dir),
        "hic_reads_dir": str(args.hic_reads_dir),
        "groups": {
            group: {
                name: file_signature(path) for name, path in paths.items()
            }
            for group, paths in group_inputs.items()
        },
    }
    return {
        "runtime": {
            "temp_dir": str(args.temp_dir) if args.temp_dir is not None else None,
        },
        "stages": {
            "scaffold": scaffold_stage_configuration(args, inputs),
        },
    }


def configuration(args):
    def signatures(names):
        return {
            name: file_signature(getattr(args, name))
            for name in names
            if getattr(args, name) is not None
        }
    return {
        "runtime": {
            "temp_dir": str(args.temp_dir) if args.temp_dir is not None else None,
        },
        "stages": {
            "assign": {
                "inputs": signatures(
                    tuple(f"bam_{kind}" for kind in args.data_types)
                    + ("contig_type", "group")
                ),
                "parameters": {
                    name: getattr(args, name)
                    for name in (
                        "seed", "collapsed_policy", "unknown_contig_type_policy",
                        "groups", "chromosomes",
                        "min_group_margin", "balanced_score_tolerance",
                        "hifi_min_alignment_length", "hifi_min_identity",
                        "ont_min_alignment_length", "ont_min_identity",
                        "hic_min_alignment_length", "hic_min_identity",
                        "hic_assignment_backend",
                    )
                    if name != "hic_assignment_backend" or "hic" in args.data_types
                } | {"data_types": args.data_types},
            },
            "extract": {
                "inputs": signatures(
                    tuple(
                        name for kind in args.data_types
                        for name in ({"hifi": ("hifi",), "ont": ("ont",),
                                      "hic": ("hic1", "hic2")}[kind])
                    )
                ),
                "parameters": {
                    "ont_length": args.ont_length,
                    "ont_quality": args.ont_quality,
                    "data_types": args.data_types,
                },
                "runtime": {
                    "backend_requested": args.extract_backend,
                    "backend_resolved": getattr(
                        args, "_extract_backend", args.extract_backend
                    ),
                    "extract_threads": args.extract_threads,
                    "tools": getattr(args, "_extract_tools", {}),
                },
            },
            "assemble": {
                "inputs": {},
                "parameters": {
                    "hifiasm": args.hifiasm,
                },
                "runtime": {
                    "jobs": args.jobs,
                    "threads_per_job": args.threads_per_job,
                },
            },
            "scaffold": scaffold_stage_configuration(args),
        },
    }


def validate_resume_configuration(previous, current, completed_stages):
    previous_stages = previous.get("stages", {})
    current_stages = current.get("stages", {})
    for stage in completed_stages:
        if _semantic_stage_configuration(previous_stages.get(stage)) != (
            _semantic_stage_configuration(current_stages.get(stage))
        ):
            raise ValueError(
                f"resume parameters for completed stage {stage} differ from run_manifest.json"
            )


def load_manifest(path: Path):
    if not path.exists():
        return None
    with path.open() as handle:
        return json.load(handle)


def save_manifest(path: Path, config, completed_stages, mode="pipeline"):
    write_json(
        path,
        {
            "format_version": 2,
            "mode": mode,
            "configuration": config,
            "completed_stages": list(completed_stages),
        },
    )


def write_json_atomic(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.tmp.{os.getpid()}"
    write_json(temporary, value)
    os.replace(temporary, path)


def configuration_key(value):
    encoded = json.dumps(
        _semantic_stage_configuration(value), sort_keys=True, default=str
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _semantic_stage_configuration(value):
    if not isinstance(value, dict):
        return value
    return {
        key: _semantic_stage_configuration(item)
        for key, item in value.items()
        if key != "runtime"
    }


def ensure_checkpoint_configuration(target: Path, value, resume):
    target.mkdir(parents=True, exist_ok=True)
    path = target / ".checkpoint_configuration.json"
    current = {"key": configuration_key(value), "configuration": value}
    if path.exists():
        previous = load_manifest(path)
        if resume and previous.get("key") != current["key"]:
            raise ValueError(
                f"Checkpoint parameters differ for {target}; rerun without --resume "
                "or choose a new --output-dir"
            )
    write_json_atomic(path, current)


def checkpoint_complete(path: Path, required_paths):
    return path.exists() and all(
        value.exists() and (not value.is_file() or value.stat().st_size > 0)
        for value in required_paths
    )


def directory_has_nonempty_file(path: Path):
    return path.is_dir() and any(
        value.is_file() and value.stat().st_size > 0
        for value in path.rglob("*")
    )


def remove_sqlite_files(path: Path):
    for value in (path, Path(f"{path}-wal"), Path(f"{path}-shm")):
        if value.exists():
            value.unlink()


def _replace_file_from_temp(source: Path, target: Path):
    """Publish a file across filesystems while keeping the final rename atomic."""
    try:
        os.replace(source, target)
        return
    except OSError as exc:
        if exc.errno != errno.EXDEV:
            raise
    publishing = target.parent / f".{target.name}.publishing.{os.getpid()}"
    try:
        shutil.copy2(source, publishing)
        os.replace(publishing, target)
        source.unlink()
    finally:
        if publishing.exists():
            publishing.unlink()


def _replace_directory_from_temp(source: Path, target: Path):
    """Publish a directory across filesystems without exposing a partial target."""
    try:
        os.replace(source, target)
        return
    except OSError as exc:
        if exc.errno != errno.EXDEV:
            raise
    publishing = target.parent / f".{target.name}.publishing.{os.getpid()}"
    try:
        shutil.copytree(source, publishing)
        os.replace(publishing, target)
        shutil.rmtree(source)
    finally:
        shutil.rmtree(publishing, ignore_errors=True)


def replace_directory_atomically(target: Path, builder, temp_directory=None):
    target.parent.mkdir(parents=True, exist_ok=True)
    backup = target.parent / f".{target.name}.previous"
    if backup.exists() and not target.exists():
        os.replace(backup, target)
    elif backup.exists():
        shutil.rmtree(backup)
    temporary = Path(
        tempfile.mkdtemp(
            prefix=f".{target.name}.",
            dir=temp_directory if temp_directory is not None else target.parent,
        )
    )
    try:
        result = builder(temporary)
        if target.exists():
            os.replace(target, backup)
        _replace_directory_from_temp(temporary, target)
        if backup.exists():
            shutil.rmtree(backup)
        return result
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        if backup.exists() and not target.exists():
            os.replace(backup, target)
        raise


def assignment_parameters(args, kind):
    return AssignmentParameters(
        min_alignment_length=getattr(args, f"{kind}_min_alignment_length"),
        min_identity=getattr(args, f"{kind}_min_identity"),
        min_total_aligned_bp=getattr(args, f"{kind}_min_alignment_length"),
        min_group_margin=args.min_group_margin,
        balanced_score_tolerance=args.balanced_score_tolerance,
        collapsed_policy=args.collapsed_policy,
    )


def run_assignment_stage(
    args, model, output_model, dosage_summary, target: Path, stage_configuration
):
    reuse = resume_stage(args, "assign")
    ensure_checkpoint_configuration(target, stage_configuration, reuse)
    checkpoints = target / ".checkpoints"
    checkpoints.mkdir(exist_ok=True)
    specifications = (
        ("hifi", args.bam_hifi, args.seed, False),
        ("ont", args.bam_ont, args.seed + 1, False),
        ("hic", args.bam_hic, args.seed + 2, True),
    )
    specifications = tuple(
        value for value in specifications if value[0] in args.data_types
    )
    type_summaries = {}
    databases = {}
    for index, (kind, bam_path, seed, is_hic) in enumerate(specifications, 1):
        LOGGER.info(
            "[assign %d/%d] %s: checking checkpoint",
            index, len(specifications), kind.upper(),
        )
        complete_path = checkpoints / f"{kind}.complete.json"
        type_summary_path = target / f"{kind}.summary.json"
        if kind == "hic" and args.hic_assignment_backend == "fast-v1":
            assignment_directory = target / "hic.fast_v1"
            required = [
                assignment_directory / "summary.json",
                *(
                    assignment_directory / f"{group}.read_ids.txt"
                    for group in model.groups
                ),
                type_summary_path,
            ]
            if reuse and complete_path.exists() and all(
                value.exists() for value in required
            ):
                LOGGER.info(
                    "[assign %d/%d] HIC: fast-v1 checkpoint complete; skipped",
                    index, len(specifications),
                )
                type_summaries[kind] = load_manifest(type_summary_path)
                databases[kind] = assignment_directory
                continue

            parameters = assignment_parameters(args, kind)
            temp_stage = stage_temp_directory(args, "01.assignments")

            def build_fast_v1(directory):
                return assign_hic_fast_v1(
                    bam_path, directory, model, parameters, seed,
                    ProgressLog("assign:hic:fast-v1", args.progress_every),
                    progress_every=args.progress_every,
                )

            fast_summary = replace_directory_atomically(
                assignment_directory, build_fast_v1, temp_stage
            )
            type_summary = {
                "backend": "fast-v1",
                "bam": fast_summary["bam_stats"],
                "assignments": fast_summary["assignments"],
                "analysis_groups": list(model.groups),
                "output_groups": list(output_model.groups),
                "read_id_directory": str(assignment_directory),
                "semantics": fast_summary["semantics"],
            }
            write_json_atomic(type_summary_path, type_summary)
            write_json_atomic(
                complete_path,
                {"complete": True, "kind": kind, "backend": "fast-v1"},
            )
            type_summaries[kind] = type_summary
            databases[kind] = assignment_directory
            LOGGER.info(
                "[assign %d/%d] HIC fast-v1 complete: %s assigned, %s deferred",
                index, len(specifications),
                f"{fast_summary['assignments']['assigned']:,}",
                f"{fast_summary['assignments']['unassigned']:,}",
            )
            continue

        decision_database = target / f"{kind}.assignments.sqlite"
        assignment_table = target / f"{kind}.assignments.tsv.gz"
        required = (decision_database, assignment_table, type_summary_path)
        if reuse and checkpoint_complete(complete_path, required):
            LOGGER.info(
                "[assign %d/%d] %s: checkpoint complete; skipped",
                index, len(specifications), kind.upper(),
            )
            type_summaries[kind] = load_manifest(type_summary_path)
            databases[kind] = decision_database
            continue

        parameters = assignment_parameters(args, kind)
        temp_stage = stage_temp_directory(args, "01.assignments")
        evidence_database = (
            temp_stage / f"{kind}.evidence.sqlite"
            if temp_stage is not None
            else checkpoints / f"{kind}.evidence.sqlite"
        )
        if not reuse:
            remove_sqlite_files(evidence_database)
        key = configuration_key(
            {"stage": stage_configuration, "kind": kind, "bam": str(bam_path)}
        )
        evidence = collect_evidence_disk(
            bam_path, model, parameters, evidence_database, key,
            hic=is_hic,
            progress=ProgressLog(f"assign:{kind}:BAM", args.progress_every),
            progress_every=args.progress_every,
        )
        targets = group_target_lengths(model, evidence["reference_lengths"])
        assignment_progress = ProgressLog(
            f"assign:{kind}:reads", max(args.progress_every // 10, 10_000)
        )
        if is_hic:
            assignment_summary = assign_hic_pairs_disk(
                evidence_database, decision_database, model, targets,
                parameters, seed, assignment_progress,
            )
        else:
            assignment_summary = assign_long_reads_disk(
                evidence_database, decision_database, model, targets,
                parameters, seed, assignment_progress,
            )
        temporary_table = (
            temp_stage / f".{kind}.assignments.{os.getpid()}.tmp.tsv.gz"
            if temp_stage is not None
            else target / f".{kind}.assignments.tmp.tsv.gz"
        )
        write_assignments_database(temporary_table, decision_database)
        _replace_file_from_temp(temporary_table, assignment_table)
        type_summary = {
            "bam": evidence["stats"],
            "resumed_from_bam_record": evidence["resumed_from_record"],
            "assignments": assignment_summary,
            "analysis_groups": list(model.groups),
            "output_groups": list(output_model.groups),
        }
        write_json_atomic(type_summary_path, type_summary)
        write_json_atomic(
            complete_path,
            {"complete": True, "kind": kind, "configuration_key": key},
        )
        remove_sqlite_files(evidence_database)
        prune_empty_temp_workspace(args)
        type_summaries[kind] = type_summary
        databases[kind] = decision_database
        LOGGER.info(
            "[assign %d/%d] %s complete: %s assigned, %s deferred",
            index, len(specifications), kind.upper(),
            f"{assignment_summary['assigned']:,}",
            f"{assignment_summary['unassigned']:,}",
        )

    summary = {
        "group_model": {
            "analysis_groups": len(model.groups),
            "output_groups": len(output_model.groups),
            "assigned_unitigs": len(model.unitig_groups),
            "ploidy": model.ploidy,
            "dosage_validation": dosage_summary,
        },
        "parameters": {
            kind: vars(assignment_parameters(args, kind))
            for kind in args.data_types
        },
        "data_types": args.data_types,
        **type_summaries,
    }
    write_json_atomic(target / "assignment_summary.json", summary)
    return databases


def load_assignment_stage(
    target: Path, data_types=("hifi", "ont", "hic"), hic_backend="fast-v1"
):
    databases = {
        kind: (
            target / "hic.fast_v1"
            if kind == "hic" and hic_backend == "fast-v1"
            else target / f"{kind}.assignments.sqlite"
        )
        for kind in data_types
    }
    missing = [str(path) for path in databases.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            f"Missing assignment checkpoint databases: {', '.join(missing)}"
        )
    return databases


def validate_missing(summary, label):
    missing = summary.get("missing_reads", summary.get("missing_pairs", 0))
    if missing:
        raise ValueError(f"{label} FASTQ is missing {missing} assigned read IDs")


def run_extraction_stage(
    args, model, databases, target: Path, stage_configuration
):
    reuse = resume_stage(args, "extract")
    ensure_checkpoint_configuration(target, stage_configuration, reuse)
    checkpoints = target / ".checkpoints"
    checkpoints.mkdir(exist_ok=True)
    summaries = {}
    specifications = (
        ("hifi", args.hifi, None, "HiFi"),
        ("ont", args.ont, None, "ONT"),
        ("hic", args.hic1, args.hic2, "Hi-C"),
    )
    specifications = tuple(
        value for value in specifications if value[0] in args.data_types
    )
    for index, (kind, input1, input2, suffix) in enumerate(specifications, 1):
        complete_path = checkpoints / f"{kind}.complete.json"
        summary_path = target / f"{kind}.summary.json"
        if kind == "hic":
            expected = [
                target / f"{group}.Hi-C.{mate}.fq.gz"
                for group in model.groups for mate in (1, 2)
            ]
        else:
            expected = [target / f"{group}.{suffix}.fq.gz" for group in model.groups]
        if reuse and checkpoint_complete(
            complete_path, [summary_path, *expected]
        ):
            LOGGER.info(
                "[extract %d/%d] %s checkpoint complete; skipped",
                index, len(specifications), kind.upper(),
            )
            summaries[kind] = load_manifest(summary_path)
            continue
        backend = (
            args._extract_backend
            if kind in ("hifi", "ont")
            else (
                "fast-v1-seqkit"
                if args.hic_assignment_backend == "fast-v1"
                else "python-paired"
            )
        )
        LOGGER.info(
            "[extract %d/%d] %s started; backend=%s",
            index, len(specifications), kind.upper(), backend,
        )
        temp_stage = stage_temp_directory(args, "02.reads")
        temporary = Path(
            tempfile.mkdtemp(
                prefix=f".{kind}.extract.",
                dir=temp_stage if temp_stage is not None else target,
            )
        )
        try:
            progress = ProgressLog(
                f"extract:{kind}", max(args.progress_every // 10, 10_000)
            )
            if kind == "hic":
                if args.hic_assignment_backend == "fast-v1":
                    summary = dispatch_paired_fastq_fast_v1(
                        input1, input2, databases[kind], model.groups, temporary,
                        seqkit=args._extract_tools["seqkit"],
                        pigz=args._extract_tools["pigz"],
                        gawk=args._extract_tools["gawk"],
                        threads=args.extract_threads,
                        jobs=args.jobs,
                        progress=progress,
                    )
                else:
                    summary = dispatch_paired_fastq_database(
                        input1, input2, databases[kind], model.groups, temporary,
                        progress,
                    )
            elif args._extract_backend == "seqkit":
                summary = dispatch_single_fastq_seqkit(
                    input1, databases[kind], model.groups, temporary, suffix,
                    min_length=args.ont_length if kind == "ont" else 0,
                    min_mean_quality=args.ont_quality if kind == "ont" else 0,
                    seqkit=args._extract_tools["seqkit"],
                    pigz=args._extract_tools["pigz"],
                    gawk=args._extract_tools["gawk"],
                    threads=args.extract_threads,
                    progress=progress,
                    progress_every=max(args.progress_every // 10, 10_000),
                )
            else:
                summary = dispatch_single_fastq_database(
                    input1, databases[kind], model.groups, temporary, suffix,
                    min_length=args.ont_length if kind == "ont" else 0,
                    min_mean_quality=args.ont_quality if kind == "ont" else 0,
                    progress=progress,
                )
            validate_missing(summary, suffix)
            for output in temporary.iterdir():
                _replace_file_from_temp(output, target / output.name)
            write_json_atomic(summary_path, summary)
            write_json_atomic(complete_path, {"complete": True, "kind": kind})
            summaries[kind] = summary
        finally:
            shutil.rmtree(temporary, ignore_errors=True)
            prune_empty_temp_workspace(args)
        LOGGER.info(
            "[extract %d/%d] %s complete",
            index, len(specifications), kind.upper(),
        )
    write_json_atomic(target / "extraction_summary.json", summaries)
    return summaries


def command_path(command):
    path = Path(command)
    if path.parent != Path("."):
        if not path.exists():
            raise FileNotFoundError(f"Tool does not exist: {command}")
        return str(path.resolve())
    resolved = shutil.which(command)
    if resolved is None:
        raise FileNotFoundError(f"Required tool is not on PATH: {command}")
    return resolved


def configure_extraction_backend(args):
    """Resolve optional fast tools once and record the actual backend."""
    args._extract_backend = "python"
    args._extract_tools = {}
    fast_hic_requires_tools = (
        "hic" in args.data_types
        and args.hic_assignment_backend == "fast-v1"
        and STAGES.index(args.stop_after) >= STAGES.index("extract")
    )
    if args.extract_backend == "python" and not fast_hic_requires_tools:
        LOGGER.info("FASTQ extraction backend: python (requested)")
        return
    requested = {
        "seqkit": args.seqkit,
        "pigz": args.pigz,
        "gawk": args.gawk,
    }
    resolved = {}
    missing = []
    for name, command in requested.items():
        try:
            resolved[name] = command_path(command)
        except FileNotFoundError as exc:
            missing.append(str(exc))
    if missing:
        if args.extract_backend == "seqkit" or fast_hic_requires_tools:
            raise FileNotFoundError("; ".join(missing))
        LOGGER.warning(
            "FASTQ extraction backend: python; fast tools unavailable: %s",
            "; ".join(missing),
        )
        return
    args._extract_backend = (
        "python" if args.extract_backend == "python" else "seqkit"
    )
    args._extract_tools = resolved
    LOGGER.info(
        "FASTQ extraction tools ready: long-read-backend=%s "
        "(seqkit=%s, gawk=%s, pigz=%s, input_threads=%d)",
        args._extract_backend,
        resolved["seqkit"], resolved["gawk"], resolved["pigz"],
        args.extract_threads,
    )


def fastq_has_records(path: Path):
    import gzip
    if not path.exists():
        return False
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return bool(handle.readline())


def run_logged(command, cwd: Path, stdout_path: Path, stderr_path: Path):
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        subprocess.run(command, cwd=cwd, stdout=stdout, stderr=stderr, check=True)


def run_parallel(items, jobs, function, label):
    items = list(items)
    total = len(items)
    if not items:
        return
    completed = 0
    with ThreadPoolExecutor(max_workers=min(jobs, max(len(items), 1))) as executor:
        futures = {executor.submit(function, item): item for item in items}
        try:
            for future in as_completed(futures):
                item = futures[future]
                future.result()
                completed += 1
                LOGGER.info(
                    "%s: %d/%d (%.1f%%) complete; latest=%s",
                    label, completed, total, 100.0 * completed / total, item,
                )
        except Exception:
            for future in futures:
                future.cancel()
            raise


def run_assembly_stage(
    args, model, reads_directory: Path, target: Path,
    stage_configuration=None,
):
    hifiasm = command_path(args.hifiasm)
    reuse = resume_stage(args, "assemble")
    ensure_checkpoint_configuration(
        target, stage_configuration or {"hifiasm": hifiasm},
        reuse,
    )
    checkpoints = target / ".checkpoints"
    checkpoints.mkdir(exist_ok=True)

    def assemble(group):
        complete = checkpoints / f"{group}.complete.json"
        final_gfa = target / group / f"{group}.asm.bp.p_ctg.gfa"
        if reuse and checkpoint_complete(complete, [final_gfa]):
            LOGGER.info("assemble checkpoint complete; skipped %s", group)
            return

        def build(group_directory):
            hifi = reads_directory / f"{group}.HiFi.fq.gz"
            ont = reads_directory / f"{group}.ONT.fq.gz"
            if not fastq_has_records(hifi):
                raise ValueError(f"No assigned HiFi reads for {group}")
            prefix = group_directory / f"{group}.asm"
            command = [hifiasm, "-t", str(args.threads_per_job), "-o", str(prefix)]
            if fastq_has_records(ont):
                command.extend(["--ul", str(ont)])
            command.append(str(hifi))
            run_logged(
                command, group_directory, group_directory / "log_hifiasm_out",
                group_directory / "log_hifiasm_err",
            )
            gfa = group_directory / f"{group}.asm.bp.p_ctg.gfa"
            if not gfa.exists() or gfa.stat().st_size == 0:
                raise RuntimeError(f"hifiasm did not create a nonempty p_ctg GFA for {group}")
        replace_directory_atomically(
            target / group, build, stage_temp_directory(args, "03.assembly")
        )
        write_json_atomic(
            complete,
            {"complete": True, "group": group, "gfa": str(final_gfa)},
        )

    groups = sorted(model.groups, key=group_sort_key)
    run_parallel(groups, args.jobs, assemble, "assemble groups")
    summary = {
        "groups": len(groups), "jobs": args.jobs,
        "threads_per_job": args.threads_per_job, "hifiasm": hifiasm,
    }
    write_json_atomic(target / "assembly_summary.json", summary)
    return summary


def gfa_to_fasta(gfa: Path, fasta: Path):
    records = 0
    with gfa.open() as source, fasta.open("w") as destination:
        for line_number, line in enumerate(source, 1):
            if not line.startswith("S\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[2] == "*":
                raise ValueError(f"GFA S record lacks sequence at {gfa}:{line_number}")
            destination.write(f">{fields[1]}\n{fields[2]}\n")
            records += 1
    if records == 0:
        raise ValueError(f"No sequence-bearing S records in {gfa}")
    return records


def run_pipe(commands, cwd: Path, stderr_prefix: Path):
    processes = []
    stderr_handles = []
    previous_stdout = None
    try:
        for index, command in enumerate(commands):
            stderr_handle = Path(f"{stderr_prefix}.{index + 1}.err").open("w")
            stderr_handles.append(stderr_handle)
            process = subprocess.Popen(
                command, cwd=cwd, stdin=previous_stdout,
                stdout=(
                    subprocess.PIPE if index < len(commands) - 1 else subprocess.DEVNULL
                ),
                stderr=stderr_handle,
            )
            if previous_stdout is not None:
                previous_stdout.close()
            previous_stdout = process.stdout
            processes.append(process)
        if previous_stdout is not None:
            previous_stdout.close()
        return_codes = [process.wait() for process in processes]
        for command, code in zip(commands, return_codes):
            if code != 0:
                raise subprocess.CalledProcessError(code, command)
    finally:
        for handle in stderr_handles:
            handle.close()


def run_scaffold_stage(
    args, model, reads_directory, assembly_directory, target,
    stage_configuration=None,
):
    tools = {
        "bwa": command_path(args.bwa), "samtools": command_path(args.samtools),
        "samblaster": command_path(args.samblaster),
        "filter_bam": command_path(args.filter_bam),
        "haphic": command_path(args.haphic),
    }
    reuse = resume_stage(args, "scaffold")
    ensure_checkpoint_configuration(target, stage_configuration or tools, reuse)
    checkpoints = target / ".checkpoints"
    checkpoints.mkdir(exist_ok=True)

    def scaffold(group):
        complete = checkpoints / f"{group}.complete.json"
        final_build = target / group / "02.haphic" / "04.build"
        if (
            reuse
            and checkpoint_complete(complete, [final_build])
            and directory_has_nonempty_file(final_build)
        ):
            LOGGER.info("scaffold checkpoint complete; skipped %s", group)
            return

        def build(group_directory):
            mapping_directory = group_directory / "01.hic_mapping"
            haphic_directory = group_directory / "02.haphic"
            mapping_directory.mkdir(parents=True)
            haphic_directory.mkdir()
            gfa = assembly_directory / group / f"{group}.asm.bp.p_ctg.gfa"
            fasta = mapping_directory / f"{group}.p_ctg.fa"
            gfa_to_fasta(gfa, fasta)
            hic1 = reads_directory / f"{group}.Hi-C.1.fq.gz"
            hic2 = reads_directory / f"{group}.Hi-C.2.fq.gz"
            if not fastq_has_records(hic1) or not fastq_has_records(hic2):
                raise ValueError(f"No assigned Hi-C pairs for {group}")
            run_logged(
                [tools["bwa"], "index", str(fasta)], mapping_directory,
                mapping_directory / "log_bwa_index_out",
                mapping_directory / "log_bwa_index_err",
            )
            bam = mapping_directory / "HiC.bam"
            run_pipe(
                [
                    [tools["bwa"], "mem", "-5SP", "-t", str(args.threads_per_job),
                     str(fasta), str(hic1), str(hic2)],
                    [tools["samblaster"]],
                    [tools["samtools"], "view", "-@", str(args.threads_per_job),
                     "-S", "-h", "-b", "-F", "3340", "-o", str(bam), "-"],
                ], mapping_directory, mapping_directory / "hic_mapping_pipe",
            )
            filtered = mapping_directory / "HiC.filtered.bam"
            run_pipe(
                [
                    [tools["filter_bam"], str(bam), "1", "--nm", str(args.haphic_nm),
                     "--threads", str(args.threads_per_job)],
                    [tools["samtools"], "view", "-b", "-@", str(args.threads_per_job),
                     "-o", str(filtered), "-"],
                ], mapping_directory, mapping_directory / "hic_filter_pipe",
            )
            run_logged(
                [tools["haphic"], "pipeline", str(fasta), str(filtered), "1",
                 "--threads", str(args.threads_per_job),
                 "--processes", str(args.haphic_processes), "--Nx", str(args.haphic_nx)],
                haphic_directory, haphic_directory / "log_haphic_out",
                haphic_directory / "log_haphic_err",
            )
            build_directory = haphic_directory / "04.build"
            if not directory_has_nonempty_file(build_directory):
                raise RuntimeError(
                    f"HapHiC did not create nonempty outputs in 04.build for {group}"
                )
            if args.run_juicebox:
                script = build_directory / "juicebox.sh"
                if not script.exists():
                    raise RuntimeError(f"HapHiC did not create juicebox.sh for {group}")
                run_logged(
                    ["bash", str(script)], build_directory,
                    build_directory / "log_juicebox_out",
                    build_directory / "log_juicebox_err",
                )
        replace_directory_atomically(
            target / group, build, stage_temp_directory(args, "04.scaffold")
        )
        write_json_atomic(
            complete,
            {"complete": True, "group": group, "build": str(final_build)},
        )

    groups = sorted(model.groups, key=group_sort_key)
    run_parallel(groups, args.jobs, scaffold, "scaffold groups")
    summary = {
        "groups": len(groups), "jobs": args.jobs,
        "threads_per_job": args.threads_per_job,
        "haphic_processes": args.haphic_processes, "tools": tools,
    }
    write_json_atomic(target / "scaffold_summary.json", summary)
    return summary


def run_scaffold_only(args, output: Path):
    manifest_path = output / "run_manifest.json"
    previous = load_manifest(manifest_path)
    if previous is not None and previous.get("mode", "pipeline") != "scaffold-only":
        raise ValueError(
            "--scaffold-only output contains a non-scaffold run_manifest.json; "
            "choose a new --output-dir"
        )

    full_model = parse_group_file(args.group)
    _, output_model = select_group_models(full_model, args)
    group_inputs = scaffold_only_group_inputs(args, output_model.groups)
    config = scaffold_only_configuration(args, group_inputs)
    completed = (
        list(previous.get("completed_stages", []))
        if args.resume and previous is not None else []
    )
    if args.rerun_from == "scaffold":
        completed = []
    if args.resume and previous is not None:
        validate_resume_configuration(
            previous.get("configuration", {}), config, completed
        )

    chromosomes = sorted(set(output_model.group_chromosome.values()))
    LOGGER.info(
        "scaffold-only selection: %d groups, chromosomes=%s",
        len(output_model.groups), ",".join(chromosomes),
    )
    LOGGER.info("scaffold-only assembly input: %s", args.assembly_dir)
    LOGGER.info("scaffold-only Hi-C input: %s", args.hic_reads_dir)
    save_manifest(manifest_path, config, completed, mode="scaffold-only")

    scaffold_directory = output / "04.scaffold"
    LOGGER.info("[scaffold-only] scaffold started/checking group checkpoints")
    run_scaffold_stage(
        args, output_model, args.hic_reads_dir, args.assembly_dir,
        scaffold_directory, config["stages"]["scaffold"],
    )
    completed = ["scaffold"]
    save_manifest(manifest_path, config, completed, mode="scaffold-only")
    LOGGER.info("[scaffold-only] scaffold complete")
    prune_empty_temp_workspace(args)
    return {"completed_stages": completed}


def run(args):
    resolved_inputs(args)
    validate_args(args)
    configure_temp_directory(args)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    log_path = output / "phase_reads.log"
    if not any(
        isinstance(handler, logging.FileHandler)
        and Path(handler.baseFilename) == log_path
        for handler in LOGGER.handlers
    ):
        file_handler = logging.FileHandler(log_path, mode="a")
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        )
        LOGGER.addHandler(file_handler)
    LOGGER.info(
        "phase_reads v2 started; mode=%s stop_after=%s resume=%s rerun_from=%s",
        "scaffold-only" if args.scaffold_only else "pipeline",
        args.stop_after, args.resume, args.rerun_from or "none",
    )
    if args.temp_dir is not None:
        LOGGER.info("temporary workspace: %s", args._temp_workspace)
    if args.scaffold_only:
        return run_scaffold_only(args, output)
    configure_extraction_backend(args)
    manifest_path = output / "run_manifest.json"
    config = configuration(args)
    previous = load_manifest(manifest_path)
    completed = list(previous.get("completed_stages", [])) if args.resume and previous else []
    if args.rerun_from:
        completed = [
            stage for stage in completed
            if STAGES.index(stage) < STAGES.index(args.rerun_from)
        ]
    if args.resume and previous is not None:
        validate_resume_configuration(
            previous.get("configuration", {}), config, completed
        )
    full_model = parse_group_file(args.group)
    dosage_summary = validate_group_dosage(
        full_model, parse_contig_types(args.contig_type),
        args.unknown_contig_type_policy,
    )
    model, output_model = select_group_models(full_model, args)
    LOGGER.info(
        "selection: %d analysis groups, %d output groups, chromosomes=%s",
        len(model.groups), len(output_model.groups),
        ",".join(sorted(set(model.group_chromosome.values()))),
    )
    assignment_directory = output / "01.assignments"
    assignment_required = [assignment_directory / "assignment_summary.json"]
    for kind in args.data_types:
        if kind == "hic" and args.hic_assignment_backend == "fast-v1":
            fast_directory = assignment_directory / "hic.fast_v1"
            assignment_required.extend(
                [
                    assignment_directory / "hic.summary.json",
                    fast_directory / "summary.json",
                    *(
                        fast_directory / f"{group}.read_ids.txt"
                        for group in model.groups
                    ),
                ]
            )
        else:
            assignment_required.extend(
                assignment_directory / f"{kind}.{suffix}"
                for suffix in ("assignments.tsv.gz", "assignments.sqlite")
            )
    assignment_files_exist = all(path.exists() for path in assignment_required)
    if args.resume and "assign" in completed and assignment_files_exist:
        decisions = load_assignment_stage(
            assignment_directory, args.data_types, args.hic_assignment_backend
        )
        LOGGER.info("[stage 1/4] assign checkpoint complete; loading disk indexes")
    else:
        LOGGER.info("[stage 1/4] assign started")
        decisions = run_assignment_stage(
            args, model, output_model, dosage_summary, assignment_directory,
            config["stages"]["assign"],
        )
        completed = ["assign"]
        save_manifest(manifest_path, config, completed)
        LOGGER.info("[stage 1/4] assign complete")
    if args.stop_after == "assign":
        prune_empty_temp_workspace(args)
        return {"completed_stages": completed}
    reads_directory = output / "02.reads"
    LOGGER.info("[stage 2/4] extract started/checking checkpoints")
    run_extraction_stage(
        args, output_model, decisions, reads_directory,
        config["stages"]["extract"],
    )
    completed = ["assign", "extract"]
    save_manifest(manifest_path, config, completed)
    LOGGER.info("[stage 2/4] extract complete")
    if args.stop_after == "extract":
        prune_empty_temp_workspace(args)
        return {"completed_stages": completed}
    assembly_directory = output / "03.assembly"
    LOGGER.info("[stage 3/4] assemble started/checking group checkpoints")
    run_assembly_stage(
        args, output_model, reads_directory, assembly_directory,
        config["stages"]["assemble"],
    )
    completed = ["assign", "extract", "assemble"]
    save_manifest(manifest_path, config, completed)
    LOGGER.info("[stage 3/4] assemble complete")
    if args.stop_after == "assemble":
        prune_empty_temp_workspace(args)
        return {"completed_stages": completed}
    scaffold_directory = output / "04.scaffold"
    LOGGER.info("[stage 4/4] scaffold started/checking group checkpoints")
    run_scaffold_stage(
        args, output_model, reads_directory, assembly_directory,
        scaffold_directory, config["stages"]["scaffold"],
    )
    completed = list(STAGES)
    save_manifest(manifest_path, config, completed)
    LOGGER.info("[stage 4/4] scaffold complete")
    prune_empty_temp_workspace(args)
    return {"completed_stages": completed}


def main(argv=None):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    args = parse_args(argv)
    summary = run(args)
    print(f"phase_reads completed stages: {','.join(summary['completed_stages'])}")


if __name__ == "__main__":
    main()
