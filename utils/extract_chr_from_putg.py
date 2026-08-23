#!/usr/bin/env python3
"""Assign selected unitigs to chromosomes and stream FASTA records.

The input PAF is expected to be the output of find_collinear_chains.py. That
step has already selected one reference target per accepted unitig,
so this script validates that decision instead of recomputing it from fragile
PAF tag positions or cumulative CIGAR lengths.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Optional


@dataclass
class ChainEvidence:
    query_length: int
    target: str
    target_length: int
    strand: str
    query_intervals: list[tuple[int, int]] = field(default_factory=list)
    matches: int = 0
    block_length: int = 0
    records: int = 0
    placement_modes: set[str] = field(default_factory=set)


@dataclass(frozen=True)
class Assignment:
    unitig: str
    query_length: Optional[int]
    target: Optional[str]
    strand: Optional[str]
    status: str
    reason: str
    placement_mode: Optional[str] = None
    records: int = 0
    query_aligned_bp: int = 0
    query_span_bp: int = 0
    query_coverage: float = 0.0
    query_span_coverage: float = 0.0
    matches: int = 0
    block_length: int = 0
    identity: float = 0.0
    reference_margin: float = 0.0


def parse_optional_tags(fields: list[str]) -> dict[str, tuple[str, str]]:
    tags = {}
    for field_value in fields:
        parts = field_value.split(":", 2)
        if len(parts) != 3 or len(parts[0]) != 2:
            continue
        tag, value_type, value = parts
        previous = tags.get(tag)
        if previous is not None and previous != (value_type, value):
            raise ValueError(f"Conflicting values for PAF tag {tag}")
        tags[tag] = (value_type, value)
    return tags


def union_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current_id = None
    current_length = 0
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_length
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                current_id = header.split()[0]
                if current_id in lengths:
                    raise ValueError(f"Duplicate FASTA ID {current_id} in {path}")
                current_length = 0
            elif line.strip():
                if current_id is None:
                    raise ValueError(f"Sequence before first FASTA header at {path}:{line_number}")
                current_length += len(line.strip())
    if current_id is not None:
        if current_id in lengths:
            raise ValueError(f"Duplicate FASTA ID {current_id} in {path}")
        lengths[current_id] = current_length
    if not lengths:
        raise ValueError(f"No FASTA records found in {path}")
    return lengths


def validate_output_name(target: str) -> None:
    if not target or target in {".", ".."} or Path(target).name != target:
        raise ValueError(f"Reference target cannot be used as an output filename: {target!r}")


def parse_collinear_paf(path: Path):
    candidates: dict[str, dict[str, ChainEvidence]] = defaultdict(dict)
    all_query_lengths: dict[str, int] = {}
    secondary_queries: set[str] = set()
    record_count = 0
    primary_record_count = 0

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            record_count += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF record at {path}:{line_number}")
            try:
                query = fields[0]
                query_length = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                strand = fields[4]
                target = fields[5]
                target_length = int(fields[6])
                target_start = int(fields[7])
                target_end = int(fields[8])
                matches = int(fields[9])
                block_length = int(fields[10])
            except ValueError as exc:
                raise ValueError(f"Invalid numeric PAF field at {path}:{line_number}") from exc

            if not (0 <= query_start < query_end <= query_length):
                raise ValueError(f"Invalid query coordinates at {path}:{line_number}")
            if not (0 <= target_start < target_end <= target_length):
                raise ValueError(f"Invalid target coordinates at {path}:{line_number}")
            if strand not in {"+", "-"}:
                raise ValueError(f"Invalid PAF strand at {path}:{line_number}: {strand}")
            previous_length = all_query_lengths.setdefault(query, query_length)
            if previous_length != query_length:
                raise ValueError(f"Inconsistent query length for {query} in {path}")

            tags = parse_optional_tags(fields[12:])
            tp = tags.get("tp")
            if tp is not None and tp == ("A", "S"):
                secondary_queries.add(query)
                continue

            primary_record_count += 1
            key = target
            evidence = candidates[query].get(key)
            if evidence is None:
                evidence = ChainEvidence(query_length, target, target_length, strand)
                candidates[query][key] = evidence
            elif evidence.target_length != target_length:
                raise ValueError(f"Inconsistent target length for {target} in {path}")
            elif evidence.strand != strand:
                evidence.strand = "mixed"
            evidence.query_intervals.append((query_start, query_end))
            evidence.matches += matches
            evidence.block_length += block_length
            evidence.records += 1
            placement = tags.get("pv")
            if placement is not None:
                evidence.placement_modes.add(placement[1])

    return candidates, all_query_lengths, secondary_queries, {
        "records": record_count,
        "primary_records": primary_record_count,
        "secondary_records": record_count - primary_record_count,
    }


def select_chromosomes(
    reference: Optional[Path],
    candidates: dict[str, dict[tuple[str, str], ChainEvidence]],
    chromosome_count: int,
):
    if chromosome_count <= 0:
        raise ValueError("--chr_num must be greater than zero")
    if reference is not None:
        lengths = read_fasta_lengths(reference)
        source = "reference_fasta"
        for query_candidates in candidates.values():
            for evidence in query_candidates.values():
                reference_length = lengths.get(evidence.target)
                if reference_length is None:
                    raise ValueError(
                        f"PAF target {evidence.target} is absent from the reference FASTA"
                    )
                if reference_length != evidence.target_length:
                    raise ValueError(
                        f"PAF target length for {evidence.target} ({evidence.target_length}) "
                        f"does not match the reference FASTA ({reference_length})"
                    )
    else:
        lengths = {}
        for query_candidates in candidates.values():
            for evidence in query_candidates.values():
                previous = lengths.setdefault(evidence.target, evidence.target_length)
                if previous != evidence.target_length:
                    raise ValueError(f"Inconsistent target length for {evidence.target}")
        source = "paf_targets"
    if not lengths:
        raise ValueError("No reference targets are available for chromosome selection")
    selected = dict(
        sorted(lengths.items(), key=lambda item: (-item[1], item[0]))[:chromosome_count]
    )
    for target in selected:
        validate_output_name(target)
    return selected, source


def build_assignments(
    candidates: dict[str, dict[tuple[str, str], ChainEvidence]],
    all_query_lengths: dict[str, int],
    secondary_queries: set[str],
    selected_chromosomes: set[str],
):
    assignments: dict[str, Assignment] = {}
    for query, query_length in all_query_lengths.items():
        query_candidates = candidates.get(query, {})
        if not query_candidates:
            reason = "secondary_only" if query in secondary_queries else "no_primary_records"
            assignments[query] = Assignment(query, query_length, None, None, "unassigned", reason)
            continue
        if len(query_candidates) != 1:
            assignments[query] = Assignment(
                query, query_length, None, None, "unassigned", "ambiguous_target"
            )
            continue

        evidence = next(iter(query_candidates.values()))
        merged = union_intervals(evidence.query_intervals)
        aligned_bp = sum(end - start for start, end in merged)
        span_bp = merged[-1][1] - merged[0][0]
        placement_mode = None
        if len(evidence.placement_modes) == 1:
            placement_mode = next(iter(evidence.placement_modes))
        elif len(evidence.placement_modes) > 1:
            assignments[query] = Assignment(
                query,
                query_length,
                evidence.target,
                evidence.strand,
                "unassigned",
                "mixed_placement_modes",
                records=evidence.records,
                query_aligned_bp=aligned_bp,
                query_span_bp=span_bp,
            )
            continue

        common = dict(
            placement_mode=placement_mode,
            records=evidence.records,
            query_aligned_bp=aligned_bp,
            query_span_bp=span_bp,
            query_coverage=aligned_bp / query_length,
            query_span_coverage=span_bp / query_length,
            matches=evidence.matches,
            block_length=evidence.block_length,
            identity=(evidence.matches / evidence.block_length if evidence.block_length else 0.0),
        )
        if evidence.target not in selected_chromosomes:
            assignments[query] = Assignment(
                query,
                query_length,
                evidence.target,
                evidence.strand,
                "unassigned",
                "target_not_selected_as_chromosome",
                **common,
            )
        else:
            assignments[query] = Assignment(
                query,
                query_length,
                evidence.target,
                evidence.strand,
                "assigned",
                "accepted_reference_placement",
                **common,
            )
    return assignments


def parse_chain_qc(path: Path):
    rows = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"query", "query_length", "best_target", "best_strand", "status", "reason"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Chain QC is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            query = row["query"]
            if query in rows:
                raise ValueError(f"Duplicate query {query} in chain QC {path}")
            rows[query] = row
    return rows


def optional_int(row, field_name, default=0):
    value = row.get(field_name, "")
    return int(value) if value not in {None, ""} else default


def optional_float(row, field_name, default=0.0):
    value = row.get(field_name, "")
    return float(value) if value not in {None, ""} else default


def enrich_assignments_from_chain_qc(
    assignments, qc_rows, selected_chromosomes, min_chromosome_margin
):
    for query, row in qc_rows.items():
        query_length = int(row["query_length"])
        existing = assignments.get(query)
        if row["status"] == "accepted":
            if existing is None:
                raise ValueError(f"Accepted chain {query} is absent from the collinear PAF")
            if existing.query_length != query_length:
                raise ValueError(f"Query length mismatch for {query} between PAF and chain QC")
            if existing.target != row["best_target"] or existing.strand != row["best_strand"]:
                raise ValueError(f"Selected target/strand mismatch for {query} between PAF and chain QC")
            assignments[query] = replace(
                existing,
                reference_margin=optional_float(row, "reference_margin"),
            )
            continue
        if existing is not None:
            raise ValueError(f"Rejected chain {query} is unexpectedly present in the collinear PAF")
        chain_reason = row["reason"] or "unspecified"
        target = row["best_target"] or None
        strand = row["best_strand"] or None
        reference_margin = optional_float(row, "reference_margin")
        status = "unassigned"
        reason = f"chain_rejected:{chain_reason}"
        if target is not None and target not in selected_chromosomes:
            reason = f"target_not_selected_as_chromosome:{chain_reason}"
        elif target is not None and strand is not None and reference_margin >= min_chromosome_margin:
            status = "assigned"
            reason = f"rescued_best_chromosome:{chain_reason}"
        elif target is not None:
            reason = f"chromosome_margin_below_threshold:{chain_reason}"
        assignments[query] = Assignment(
            query,
            query_length,
            target,
            strand,
            status,
            reason,
            placement_mode=row.get("placement_mode") or None,
            records=optional_int(row, "chain_alignments"),
            query_aligned_bp=optional_int(row, "query_covered_bp"),
            query_span_bp=optional_int(row, "query_span_bp"),
            query_coverage=optional_float(row, "query_coverage"),
            query_span_coverage=optional_float(row, "query_span_coverage"),
            identity=optional_float(row, "chain_identity"),
            reference_margin=reference_margin,
        )
    return assignments


ASSIGNMENT_FIELDS = [
    "unitig",
    "fasta_length",
    "paf_query_length",
    "target",
    "strand",
    "status",
    "reason",
    "placement_mode",
    "paf_records",
    "query_aligned_bp",
    "query_span_bp",
    "query_coverage",
    "query_span_coverage",
    "identity",
    "reference_margin",
]


def assignment_row(assignment: Assignment, fasta_length: int):
    return {
        "unitig": assignment.unitig,
        "fasta_length": fasta_length,
        "paf_query_length": assignment.query_length if assignment.query_length is not None else "",
        "target": assignment.target or "",
        "strand": assignment.strand or "",
        "status": assignment.status,
        "reason": assignment.reason,
        "placement_mode": assignment.placement_mode or "",
        "paf_records": assignment.records,
        "query_aligned_bp": assignment.query_aligned_bp,
        "query_span_bp": assignment.query_span_bp,
        "query_coverage": f"{assignment.query_coverage:.6f}",
        "query_span_coverage": f"{assignment.query_span_coverage:.6f}",
        "identity": f"{assignment.identity:.6f}",
        "reference_margin": f"{assignment.reference_margin:.6f}",
    }


def stream_fasta_outputs(
    fasta: Path,
    temporary_directory: Path,
    assignments: dict[str, Assignment],
):
    chromosome_handles = {}
    chromosome_stats = defaultdict(lambda: {"unitigs": 0, "bp": 0})
    unassigned_stats = {"unitigs": 0, "bp": 0}
    rows = []
    seen = set()
    current_id = None
    current_assignment = None
    current_handle = None
    current_length = 0
    unassigned_path = temporary_directory / "un_chr.fa"

    def finish_record():
        if current_id is None:
            return
        if (
            current_assignment.query_length is not None
            and current_assignment.query_length != current_length
        ):
            raise ValueError(
                f"PAF query length for {current_id} ({current_assignment.query_length}) "
                f"does not match the p_utg FASTA ({current_length})"
            )
        if current_assignment.status == "assigned":
            stats = chromosome_stats[current_assignment.target]
        else:
            stats = unassigned_stats
        stats["unitigs"] += 1
        stats["bp"] += current_length
        rows.append(assignment_row(current_assignment, current_length))

    with fasta.open() as source, unassigned_path.open("w") as unassigned_handle:
        try:
            for line_number, line in enumerate(source, 1):
                if line.startswith(">"):
                    finish_record()
                    header = line[1:].strip()
                    if not header:
                        raise ValueError(f"Empty FASTA header at {fasta}:{line_number}")
                    current_id = header.split()[0]
                    if current_id in seen:
                        raise ValueError(f"Duplicate FASTA ID {current_id} in {fasta}")
                    seen.add(current_id)
                    current_assignment = assignments.get(current_id)
                    if current_assignment is None:
                        current_assignment = Assignment(
                            current_id,
                            None,
                            None,
                            None,
                            "unassigned",
                            "no_accepted_collinear_chain",
                        )

                    if current_assignment.status == "assigned":
                        target = current_assignment.target
                        current_handle = chromosome_handles.get(target)
                        if current_handle is None:
                            path = temporary_directory / f"{target}.putg.fa"
                            current_handle = path.open("w")
                            chromosome_handles[target] = current_handle
                    else:
                        current_handle = unassigned_handle
                    current_handle.write(line if line.endswith("\n") else line + "\n")
                    current_length = 0
                elif line.strip():
                    if current_id is None:
                        raise ValueError(
                            f"Sequence before first FASTA header at {fasta}:{line_number}"
                        )
                    sequence = line.strip()
                    current_handle.write(sequence + "\n")
                    current_length += len(sequence)
            finish_record()
        finally:
            for handle in chromosome_handles.values():
                handle.close()

    missing = sorted(set(assignments) - seen)
    if missing:
        examples = ", ".join(missing[:5])
        raise ValueError(
            f"{len(missing)} PAF query IDs are missing from the p_utg FASTA; examples: {examples}"
        )
    return rows, dict(chromosome_stats), unassigned_stats


def write_metadata(
    temporary_directory: Path,
    rows,
    selected_chromosomes,
    selection_source,
    chromosome_stats,
    unassigned_stats,
    paf_stats,
    paths,
    min_chromosome_margin,
):
    assignment_path = temporary_directory / "chromosome_assignments.tsv"
    with assignment_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=ASSIGNMENT_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    with (temporary_directory / "putg_vs_mT2T.best.match.txt").open("w") as handle:
        handle.write("#contig_id\tscaffold_id\n")
        for row in rows:
            if row["status"] == "assigned":
                handle.write(f'{row["unitig"]}\t{row["target"]}\n')

    reason_counts = Counter(row["reason"] for row in rows)
    assigned_unitigs = sum(value["unitigs"] for value in chromosome_stats.values())
    assigned_bp = sum(value["bp"] for value in chromosome_stats.values())
    summary = {
        "inputs": {key: str(value.resolve()) for key, value in paths.items() if value is not None},
        "chromosome_selection_source": selection_source,
        "selected_chromosomes": selected_chromosomes,
        "paf": paf_stats,
        "parameters": {"min_chromosome_margin": min_chromosome_margin},
        "fasta": {
            "unitigs": assigned_unitigs + unassigned_stats["unitigs"],
            "bp": assigned_bp + unassigned_stats["bp"],
        },
        "assigned": {"unitigs": assigned_unitigs, "bp": assigned_bp},
        "unassigned": unassigned_stats,
        "assignment_reason_counts": dict(sorted(reason_counts.items())),
        "chromosomes": dict(sorted(chromosome_stats.items())),
    }
    with (temporary_directory / "chromosome_assignment.summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def install_outputs(temporary_directory: Path, output_directory: Path):
    generated = {path.name for path in temporary_directory.iterdir() if path.is_file()}
    for old_path in output_directory.glob("*.putg.fa"):
        if old_path.name not in generated:
            old_path.unlink()
    for obsolete_name in (
        "contig_match_ratios.csv",
        "contig_match_ratios.xlsx",
        "putg_vs_mT2T.match_above_cutoff.txt",
    ):
        obsolete = output_directory / obsolete_name
        if obsolete.exists():
            obsolete.unlink()
    for temporary_path in sorted(temporary_directory.iterdir()):
        if temporary_path.is_file():
            os.replace(temporary_path, output_directory / temporary_path.name)


def run(args):
    if not 0.0 <= args.min_chromosome_margin <= 1.0:
        raise ValueError("--min-chromosome-margin must be between zero and one")
    output_directory = args.wd.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    candidates, query_lengths, secondary_queries, paf_stats = parse_collinear_paf(args.paf)
    selected, selection_source = select_chromosomes(args.mT2T, candidates, args.chr_num)
    assignments = build_assignments(candidates, query_lengths, secondary_queries, set(selected))
    if args.chain_qc is not None:
        assignments = enrich_assignments_from_chain_qc(
            assignments,
            parse_chain_qc(args.chain_qc),
            set(selected),
            args.min_chromosome_margin,
        )

    temporary_directory = Path(tempfile.mkdtemp(prefix=".chr_seq.", dir=output_directory))
    try:
        rows, chromosome_stats, unassigned_stats = stream_fasta_outputs(
            args.p_utg, temporary_directory, assignments
        )
        summary = write_metadata(
            temporary_directory,
            rows,
            selected,
            selection_source,
            chromosome_stats,
            unassigned_stats,
            paf_stats,
            {
                "p_utg": args.p_utg,
                "paf": args.paf,
                "mT2T": args.mT2T,
                "chain_qc": args.chain_qc,
            },
            args.min_chromosome_margin,
        )
        install_outputs(temporary_directory, output_directory)
    finally:
        shutil.rmtree(temporary_directory, ignore_errors=True)
    print(
        "Chromosome extraction: "
        f'{summary["assigned"]["unitigs"]} assigned, '
        f'{summary["unassigned"]["unitigs"]} unassigned'
    )
    return summary


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Assign accepted collinear p_utg chains to mT2T chromosomes"
    )
    parser.add_argument("--wd", required=True, type=Path, help="Output directory")
    parser.add_argument("--paf", required=True, type=Path, help="Accepted collinear PAF")
    parser.add_argument("--chain-qc", type=Path, help="Collinear-chain QC with rejection reasons")
    parser.add_argument("--chr_num", type=int, default=12, help="Number of chromosomes [12]")
    parser.add_argument(
        "--min-chromosome-margin",
        type=float,
        default=0.05,
        help="Minimum (best - second) / best chain score for chromosome rescue [0.05]",
    )
    parser.add_argument("--p_utg", required=True, type=Path, help="p_utg FASTA")
    parser.add_argument(
        "--mT2T",
        type=Path,
        help="mT2T FASTA; recommended for chromosome selection independent of PAF hits",
    )
    return parser.parse_args()


def main():
    run(parse_arguments())


if __name__ == "__main__":
    main()
