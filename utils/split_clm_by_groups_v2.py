#!/usr/bin/env python3
"""Split a Hi-C CLM once across all chromosome haplotype groups."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import tempfile
from collections import Counter
from pathlib import Path


GROUP_PATTERN = re.compile(r"^(.+)_group([0-9]+)$")


def load_groups(recluster_directory: Path, ploidy: int):
    unitig_groups = {}
    unitig_chromosomes = {}
    group_paths = {}
    group_sizes = Counter()
    cluster_files = sorted(recluster_directory.glob("chr*/group.reassignment.cluster.txt"))
    if not cluster_files:
        raise ValueError(f"No chromosome reassignment cluster files in {recluster_directory}")
    for path in cluster_files:
        with path.open() as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.split()
                match = GROUP_PATTERN.match(fields[0])
                if not match:
                    raise ValueError(f"Invalid group ID at {path}:{line_number}: {fields[0]}")
                chromosome, group_text = match.groups()
                group = int(group_text)
                if group < 1 or group > ploidy:
                    raise ValueError(f"Out-of-range group at {path}:{line_number}")
                key = (chromosome, group)
                if key in group_paths:
                    raise ValueError(f"Duplicate group definition: {fields[0]}")
                group_paths[key] = path.parent / "split_clms" / f"group{group}.clm"
                for unitig in fields[1:]:
                    previous = unitig_chromosomes.setdefault(unitig, chromosome)
                    if previous != chromosome:
                        raise ValueError(
                            f"Unitig {unitig} occurs in both {previous} and {chromosome}"
                        )
                    unitig_groups.setdefault(unitig, set()).add(key)
                    group_sizes[key] += 1
    expected = len(cluster_files) * ploidy
    if len(group_paths) != expected:
        raise ValueError(
            f"Expected {expected} chromosome groups, found {len(group_paths)}"
        )
    return unitig_groups, group_paths, group_sizes, cluster_files


def strip_orientation(token: str):
    return token[:-1] if token.endswith(("+", "-")) else token


def load_cluster_file(cluster_file: Path, output_directory: Path, ploidy: int):
    unitig_groups = {}
    unitig_chromosomes = {}
    group_paths = {}
    group_sizes = Counter()
    chromosome_groups = {}
    with cluster_file.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            match = GROUP_PATTERN.match(fields[0])
            if not match:
                raise ValueError(
                    f"Invalid group ID at {cluster_file}:{line_number}: {fields[0]}"
                )
            chromosome, group_text = match.groups()
            group = int(group_text)
            if group < 1 or group > ploidy:
                raise ValueError(f"Out-of-range group at {cluster_file}:{line_number}")
            key = (chromosome, group)
            if key in group_paths:
                raise ValueError(f"Duplicate group definition: {fields[0]}")
            chromosome_groups.setdefault(chromosome, set()).add(group)
            group_paths[key] = output_directory / "split_clms" / f"{fields[0]}.clm"
            for unitig in fields[1:]:
                previous = unitig_chromosomes.setdefault(unitig, chromosome)
                if previous != chromosome:
                    raise ValueError(
                        f"Unitig {unitig} occurs in both {previous} and {chromosome}"
                    )
                unitig_groups.setdefault(unitig, set()).add(key)
                group_sizes[key] += 1
    if not group_paths:
        raise ValueError(f"No groups found in {cluster_file}")
    expected_groups = set(range(1, ploidy + 1))
    for chromosome, observed in chromosome_groups.items():
        if observed != expected_groups:
            raise ValueError(
                f"Incomplete group set for {chromosome}: expected "
                f"{sorted(expected_groups)}, observed {sorted(observed)}"
            )
    return unitig_groups, group_paths, group_sizes, [cluster_file]


def split_group_mapping(
    clm_path: Path, unitig_groups, group_paths, group_sizes,
    cluster_files, output_root: Path, input_description,
):
    temporary_root = Path(
        tempfile.mkdtemp(prefix=".split_clms.", dir=output_root)
    )
    handles = {}
    output_records = Counter()
    total_records = 0
    written_records = 0
    malformed = 0
    try:
        for chromosome, group in sorted(group_paths):
            path = temporary_root / chromosome / f"group{group}.clm"
            path.parent.mkdir(parents=True, exist_ok=True)
            handles[(chromosome, group)] = path.open("w")

        with clm_path.open() as source:
            for line_number, line in enumerate(source, 1):
                if not line.strip():
                    continue
                total_records += 1
                fields = line.split(None, 2)
                if len(fields) < 2:
                    malformed += 1
                    continue
                unitig1 = strip_orientation(fields[0])
                unitig2 = strip_orientation(fields[1])
                common = unitig_groups.get(unitig1, set()) & unitig_groups.get(unitig2, set())
                for key in sorted(common):
                    handles[key].write(line)
                    output_records[key] += 1
                    written_records += 1
        if malformed:
            raise ValueError(f"CLM contains {malformed} malformed nonempty records")
    except Exception:
        for handle in handles.values():
            handle.close()
        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
    else:
        for handle in handles.values():
            handle.close()

    try:
        expected_destinations = {path.resolve() for path in group_paths.values()}
        for directory in {path.parent for path in group_paths.values()}:
            if directory.is_dir():
                for old_path in directory.glob("*.clm"):
                    if old_path.resolve() not in expected_destinations:
                        old_path.unlink()
        for chromosome, group in sorted(group_paths):
            source = temporary_root / chromosome / f"group{group}.clm"
            destination_directory = group_paths[(chromosome, group)].parent
            destination_directory.mkdir(parents=True, exist_ok=True)
            os.replace(source, group_paths[(chromosome, group)])
        summary = {
            "input": str(clm_path.resolve()),
            **input_description,
            "cluster_files": [str(path.resolve()) for path in cluster_files],
            "input_records": total_records,
            "written_group_records": written_records,
            "groups": {
                f"{chromosome}_group{group}": {
                    "unitigs": group_sizes[(chromosome, group)],
                    "clm_records": output_records[(chromosome, group)],
                    "path": str(group_paths[(chromosome, group)].resolve()),
                }
                for chromosome, group in sorted(group_paths)
            },
        }
        summary_path = output_root / "clm_split.summary.json"
        temporary_summary = output_root / ".clm_split.summary.json.tmp"
        with temporary_summary.open("w") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary_summary, summary_path)
        return summary
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)


def split_clm(clm_path: Path, recluster_directory: Path, ploidy: int):
    unitig_groups, group_paths, group_sizes, cluster_files = load_groups(
        recluster_directory, ploidy
    )
    return split_group_mapping(
        clm_path, unitig_groups, group_paths, group_sizes, cluster_files,
        recluster_directory,
        {"recluster_directory": str(recluster_directory.resolve())},
    )


def split_cluster_file_clm(
    clm_path: Path, cluster_file: Path, output_directory: Path, ploidy: int
):
    unitig_groups, group_paths, group_sizes, cluster_files = load_cluster_file(
        cluster_file, output_directory, ploidy
    )
    return split_group_mapping(
        clm_path, unitig_groups, group_paths, group_sizes, cluster_files,
        output_directory,
        {
            "clusters_file": str(cluster_file.resolve()),
            "output_directory": str(output_directory.resolve()),
        },
    )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Split one CLM into all PHap v2 chromosome haplotype groups"
    )
    parser.add_argument("--clm", required=True, type=Path)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--recluster-dir", type=Path)
    source.add_argument("--clusters-file", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--ploidy", type=int, default=4)
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.recluster_dir:
        if args.output_dir:
            raise ValueError("--output-dir is only valid with --clusters-file")
        summary = split_clm(args.clm, args.recluster_dir.resolve(), args.ploidy)
    else:
        if not args.output_dir:
            raise ValueError("--output-dir is required with --clusters-file")
        output_directory = args.output_dir.resolve()
        output_directory.mkdir(parents=True, exist_ok=True)
        summary = split_cluster_file_clm(
            args.clm, args.clusters_file.resolve(), output_directory, args.ploidy
        )
    print(
        f"Split {summary['input_records']} CLM records across "
        f"{len(summary['groups'])} chromosome groups in one pass"
    )


if __name__ == "__main__":
    main()
