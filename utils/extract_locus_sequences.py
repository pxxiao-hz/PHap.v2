#!/usr/bin/env python3
"""Partition unitigs by strict mT2T locus decisions without assigning groups."""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path
from typing import Dict, Iterable, Tuple

from phap_core.fasta import FastaRecord, fasta_lines, read_fasta


_SAFE_LOCUS_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


def parse_decisions(path: str) -> Dict[str, Tuple[str, str, str]]:
    """Return ``unitig -> (status, assigned_locus, reason)``."""

    with Path(path).open(encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        required = {"unitig_ID", "status", "assigned_locus", "reason"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                "locus decision table must contain unitig_ID, status, "
                "assigned_locus, and reason"
            )
        decisions: Dict[str, Tuple[str, str, str]] = {}
        for line_number, row in enumerate(reader, start=2):
            unitig_id = row["unitig_ID"]
            if not unitig_id or unitig_id in decisions:
                raise ValueError(
                    f"{path}:{line_number}: invalid or duplicate unitig ID"
                )
            locus = row["assigned_locus"]
            decisions[unitig_id] = (
                row["status"],
                "" if locus == "." else locus,
                row["reason"],
            )
    if not decisions:
        raise ValueError("locus decision table contains no records")
    return decisions


def partition_records(
    records: Iterable[FastaRecord],
    decisions: Dict[str, Tuple[str, str, str]],
) -> Tuple[
    Dict[str, Tuple[FastaRecord, ...]],
    Tuple[FastaRecord, ...],
    Tuple[Tuple[str, str, str, str], ...],
]:
    """Create one disjoint destination for every source FASTA record."""

    record_list = tuple(records)
    record_ids = {record.identifier for record in record_list}
    unknown_decisions = sorted(set(decisions) - record_ids)
    if unknown_decisions:
        raise ValueError(
            "locus decisions reference IDs absent from FASTA: "
            + ", ".join(unknown_decisions[:5])
        )

    by_locus: Dict[str, list[FastaRecord]] = {}
    unplaced = []
    routing = []
    for record in sorted(record_list, key=lambda row: row.identifier):
        decision = decisions.get(record.identifier)
        if decision is None:
            status, locus, reason = (
                "unassigned",
                "",
                "missing_locus_decision",
            )
        else:
            status, locus, reason = decision
        if locus:
            if _SAFE_LOCUS_ID.fullmatch(locus) is None:
                raise ValueError(f"unsafe locus ID for output filename: {locus!r}")
            by_locus.setdefault(locus, []).append(record)
            destination = f"{locus}.putg.fa"
        else:
            unplaced.append(record)
            destination = "un_chr.fa"
        routing.append((record.identifier, destination, status, reason))

    assigned_ids = {
        record.identifier
        for locus_records in by_locus.values()
        for record in locus_records
    }
    unplaced_ids = {record.identifier for record in unplaced}
    if assigned_ids & unplaced_ids or assigned_ids | unplaced_ids != record_ids:
        raise AssertionError("FASTA locus partition must be complete and disjoint")
    return (
        {
            locus: tuple(sorted(locus_records, key=lambda row: row.identifier))
            for locus, locus_records in sorted(by_locus.items())
        },
        tuple(sorted(unplaced, key=lambda row: row.identifier)),
        tuple(routing),
    )


def write_partition(
    records: Iterable[FastaRecord],
    decisions: Dict[str, Tuple[str, str, str]],
    output_directory: str,
) -> None:
    by_locus, unplaced, routing = partition_records(records, decisions)
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    manifest = output / "locus_sequence_manifest.tsv"
    _remove_manifest_owned_locus_files(output, manifest)
    for locus, locus_records in by_locus.items():
        _atomic_write_lines(output / f"{locus}.putg.fa", fasta_lines(locus_records))
    _atomic_write_lines(output / "un_chr.fa", fasta_lines(unplaced))
    _atomic_write_lines(
        output / "locus_sequence_routing.tsv",
        [
            "source_ID\tdestination\tstatus\treason",
            *("\t".join(row) for row in routing),
        ],
    )
    _atomic_write_lines(
        manifest,
        [
            "output_file\trecord_count",
            *(
                f"{locus}.putg.fa\t{len(locus_records)}"
                for locus, locus_records in by_locus.items()
            ),
            f"un_chr.fa\t{len(unplaced)}",
        ],
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p-utg", required=True, help="Source unitig FASTA")
    parser.add_argument("--locus-decisions", required=True)
    parser.add_argument("--output-directory", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    write_partition(
        read_fasta(args.p_utg),
        parse_decisions(args.locus_decisions),
        args.output_directory,
    )


def _atomic_write_lines(path: Path, lines: Iterable[str]) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as output:
        for line in lines:
            output.write(line)
            output.write("\n")
    os.replace(temporary, path)


def _remove_manifest_owned_locus_files(output: Path, manifest: Path) -> None:
    existing = {path.name for path in output.glob("*.putg.fa")}
    if not manifest.exists():
        if existing:
            raise ValueError(
                "output directory contains locus FASTA without a PHap manifest: "
                + ", ".join(sorted(existing))
            )
        return

    rows = [
        line.split("\t")
        for line in manifest.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not rows or rows[0] != ["output_file", "record_count"]:
        raise ValueError("invalid locus sequence manifest")
    owned: set[str] = set()
    for row in rows[1:]:
        if len(row) != 2:
            raise ValueError("invalid locus sequence manifest row")
        filename = row[0]
        if filename == "un_chr.fa":
            continue
        if (
            Path(filename).name != filename
            or not filename.endswith(".putg.fa")
        ):
            raise ValueError(f"unsafe locus sequence manifest path: {filename!r}")
        owned.add(filename)
    unknown = existing - owned
    if unknown:
        raise ValueError(
            "output directory contains locus FASTA not owned by the manifest: "
            + ", ".join(sorted(unknown))
        )
    for filename in sorted(owned):
        path = output / filename
        if path.exists():
            path.unlink()


if __name__ == "__main__":
    main()
