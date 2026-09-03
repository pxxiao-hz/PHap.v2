#!/usr/bin/env python3
"""Rescue chromosome-unassigned unitigs into trusted chromosome haplotype groups."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import pickle
import re
import shutil
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

from dosage import dosage_from_contig_type

GROUP_PATTERN = re.compile(r"^(.+)_group([0-9]+)$")


def parse_fasta(path: Path, motif: str, flank: Optional[int]):
    records = {}
    order = []
    current = None
    chunks = []

    def finish():
        if current is None:
            return
        sequence = "".join(chunks)
        upper = sequence.upper()
        if flank is None:
            re_sequence = upper
        else:
            flank_size = min(flank, len(sequence) // 2)
            re_sequence = (
                upper[:flank_size] + upper[-flank_size:] if flank_size else ""
            )
        records[current] = {
            "sequence": sequence,
            "length": len(sequence),
            "re_sites": re_sequence.count(motif) + 1,
        }

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                finish()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                current = header.split()[0]
                if current in records or current in order:
                    raise ValueError(f"Duplicate FASTA ID {current} at {path}:{line_number}")
                order.append(current)
                chunks = []
            elif line.strip():
                if current is None:
                    raise ValueError(
                        f"Sequence before first FASTA header at {path}:{line_number}"
                    )
                chunks.append(line.strip())
    finish()
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records, order


def scan_fasta_metadata(path: Path):
    records = {}
    order = []
    current = None
    length = 0

    def finish():
        if current is not None:
            records[current] = length

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                finish()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                current = header.split()[0]
                if current in records or current in order:
                    raise ValueError(f"Duplicate FASTA ID {current} at {path}:{line_number}")
                order.append(current)
                length = 0
            elif line.strip():
                if current is None:
                    raise ValueError(
                        f"Sequence before first FASTA header at {path}:{line_number}"
                    )
                length += len(line.strip())
    finish()
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records, order


def parse_contig_types(path: Path):
    result = {}
    with path.open() as handle:
        header = handle.readline().split()
        try:
            id_index = header.index("contig_ID")
            type_index = header.index("contig_type")
        except ValueError as exc:
            raise ValueError(
                "contig type file must contain contig_ID and contig_type columns"
            ) from exc
        for line_number, line in enumerate(handle, 2):
            fields = line.split()
            if not fields:
                continue
            if len(fields) <= max(id_index, type_index):
                raise ValueError(f"Malformed contig type row at {path}:{line_number}")
            unitig = fields[id_index]
            if unitig in result:
                raise ValueError(f"Duplicate contig type for {unitig} at {path}:{line_number}")
            result[unitig] = fields[type_index]
    return result


def candidate_dosages(records, contig_types, ploidy, policy):
    dosage = {}
    dosage_source = {}
    unsupported = {}
    for unitig in records:
        contig_type = contig_types.get(unitig)
        value = dosage_from_contig_type(contig_type)
        reason = None
        if contig_type is None:
            reason = "missing_contig_type"
        elif value is None:
            reason = f"unsupported_contig_type:{contig_type}"
        elif value > ploidy:
            reason = f"dosage_exceeds_ploidy:{value}>{ploidy}"
        if reason is None:
            dosage[unitig] = value
            dosage_source[unitig] = "provided"
        elif policy == "error":
            raise ValueError(f"Cannot determine rescue dosage for {unitig}: {reason}")
        elif policy == "haplotig":
            dosage[unitig] = 1
            dosage_source[unitig] = "default_haplotig"
        else:
            unsupported[unitig] = reason
            dosage_source[unitig] = "deferred"
    return dosage, dosage_source, unsupported


def load_recluster_state(directory: Path, ploidy: int):
    inventory = {}
    assignments = {}
    group_members = {}
    group_fasta_paths = {}
    unassigned_fasta_paths = []
    chromosomes = []
    table_paths = sorted(directory.glob("chr*/recluster_assignments.tsv"))
    if not table_paths:
        raise ValueError(f"No chromosome reassignment tables in {directory}")

    for table_path in table_paths:
        chromosome = table_path.parent.name
        chromosomes.append(chromosome)
        with table_path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "unitig", "length", "re_sites", "contig_type", "dosage",
                "dosage_source", "status", "groups", "assignment_basis",
            }
            missing = required - set(reader.fieldnames or ())
            if missing:
                raise ValueError(
                    f"Missing columns in {table_path}: {', '.join(sorted(missing))}"
                )
            for row in reader:
                unitig = row["unitig"]
                if unitig in inventory:
                    raise ValueError(f"Unitig {unitig} occurs in multiple chromosomes")
                groups = tuple(
                    f"{chromosome}_group{int(group)}"
                    for group in row["groups"].split(",") if group
                )
                status = row["status"]
                if status not in {"assigned", "unassigned"}:
                    raise ValueError(f"Invalid status for {unitig} in {table_path}: {status}")
                record = {
                    "length": int(row["length"]),
                    "re_sites": int(row["re_sites"]),
                    "contig_type": row["contig_type"],
                    "dosage": int(row["dosage"]),
                    "dosage_source": row["dosage_source"],
                    "status": status,
                    "chromosome": chromosome,
                    "groups": groups,
                    "assignment_basis": row["assignment_basis"],
                    "source": (
                        "recluster_assigned" if status == "assigned"
                        else "recluster_deferred"
                    ),
                }
                if status == "assigned":
                    if len(groups) != record["dosage"]:
                        raise ValueError(
                            f"Recluster dosage mismatch for {unitig}: "
                            f"expected {record['dosage']}, observed {len(groups)}"
                        )
                    assignments[unitig] = groups
                elif groups:
                    raise ValueError(f"Unassigned recluster unitig {unitig} has groups")
                inventory[unitig] = record

        cluster_path = table_path.parent / "group.reassignment.cluster.txt"
        observed_memberships = defaultdict(set)
        seen_groups = set()
        with cluster_path.open() as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.split()
                match = GROUP_PATTERN.match(fields[0])
                if not match or match.group(1) != chromosome:
                    raise ValueError(
                        f"Invalid group ID at {cluster_path}:{line_number}: {fields[0]}"
                    )
                group_number = int(match.group(2))
                if group_number < 1 or group_number > ploidy or fields[0] in seen_groups:
                    raise ValueError(f"Duplicate or out-of-range group at {cluster_path}:{line_number}")
                seen_groups.add(fields[0])
                if len(fields[1:]) != len(set(fields[1:])):
                    raise ValueError(f"Duplicate member at {cluster_path}:{line_number}")
                group_members[fields[0]] = list(fields[1:])
                group_fasta_paths[fields[0]] = (
                    table_path.parent / f"group{group_number}.reassignment.fa"
                )
                for unitig in fields[1:]:
                    observed_memberships[unitig].add(fields[0])
        expected_groups = {
            f"{chromosome}_group{group}" for group in range(1, ploidy + 1)
        }
        if seen_groups != expected_groups:
            raise ValueError(f"Incomplete group set in {cluster_path}")
        expected_memberships = {
            unitig: set(groups)
            for unitig, groups in assignments.items()
            if inventory[unitig]["chromosome"] == chromosome
        }
        if dict(observed_memberships) != expected_memberships:
            raise ValueError(f"Cluster/table membership mismatch for {chromosome}")
        unassigned_path = table_path.parent / "unassigned_unitigs.fa"
        if not unassigned_path.is_file():
            raise ValueError(f"Missing recluster unassigned FASTA: {unassigned_path}")
        unassigned_fasta_paths.append((chromosome, unassigned_path))

    return {
        "inventory": inventory,
        "assignments": assignments,
        "group_members": group_members,
        "group_fasta_paths": group_fasta_paths,
        "unassigned_fasta_paths": unassigned_fasta_paths,
        "chromosomes": sorted(chromosomes),
        "table_paths": table_paths,
    }


def load_relevant_links(
    path: Path, candidates, assignments, dosage, normalization="dosage"
):
    if normalization not in {"dosage", "raw"}:
        raise ValueError(
            f"Unsupported Hi-C link normalization: {normalization}; expected dosage or raw"
        )
    with path.open("rb") as handle:
        raw_links = pickle.load(handle)
    if not isinstance(raw_links, dict):
        raise ValueError("Hi-C link pickle must contain a dictionary")
    eligible = set(assignments) | set(dosage)
    candidate_set = set(candidates)
    adjusted = defaultdict(float)
    invalid = 0
    for key, value in raw_links.items():
        if not isinstance(key, tuple) or len(key) != 2:
            invalid += 1
            continue
        unitig1, unitig2 = key
        if (
            unitig1 == unitig2
            or unitig1 not in eligible
            or unitig2 not in eligible
            or not ({unitig1, unitig2} & candidate_set)
        ):
            continue
        try:
            count = float(value)
        except (TypeError, ValueError):
            invalid += 1
            continue
        if not math.isfinite(count) or count < 0:
            invalid += 1
            continue
        pair = tuple(sorted((unitig1, unitig2)))
        if normalization == "dosage":
            count /= dosage[unitig1] * dosage[unitig2]
        adjusted[pair] += count
    if invalid:
        raise ValueError(f"Hi-C link pickle contains {invalid} invalid records")
    neighbors = defaultdict(dict)
    for (unitig1, unitig2), count in adjusted.items():
        neighbors[unitig1][unitig2] = count
        neighbors[unitig2][unitig1] = count
    return dict(neighbors), len(raw_links), len(adjusted)


def initial_group_re_sites(assignments, inventory, group_ids):
    totals = {group: 1 for group in group_ids}
    for unitig, groups in assignments.items():
        for group in groups:
            totals[group] += inventory[unitig]["re_sites"] - 1
    return totals


def score_candidate(
    unitig, assignments, neighbors, group_re_sites, chromosomes, ploidy, dosage
):
    group_links = defaultdict(float)
    for neighbor, count in neighbors.get(unitig, {}).items():
        for group in assignments.get(neighbor, ()):
            group_links[group] += count
    densities = {
        group: group_links[group] / max(group_re_sites[group], 1)
        for chromosome in chromosomes
        for group in (
            f"{chromosome}_group{number}" for number in range(1, ploidy + 1)
        )
    }
    chromosome_metrics = []
    for chromosome in chromosomes:
        local_groups = [
            f"{chromosome}_group{number}" for number in range(1, ploidy + 1)
        ]
        ranked = sorted(local_groups, key=lambda group: (-densities[group], group))
        selected = tuple(
            sorted(ranked[:dosage], key=lambda group: int(group.rsplit("group", 1)[1]))
        )
        weakest_group = min(
            selected, key=lambda group: (densities[group], group), default=None
        )
        weakest = densities[weakest_group] if weakest_group is not None else 0.0
        unselected = [group for group in local_groups if group not in selected]
        strongest_group = max(
            unselected, key=lambda group: (densities[group], group), default=None
        )
        strongest_unselected = (
            densities[strongest_group] if strongest_group is not None else 0.0
        )
        if not unselected:
            group_margin = 1.0
        elif weakest > 0:
            group_margin = (weakest - strongest_unselected) / weakest
        else:
            group_margin = 0.0
        chromosome_metrics.append({
            "chromosome": chromosome,
            "groups": selected,
            "score": sum(densities[group] for group in selected),
            "selected_links": sum(group_links[group] for group in selected),
            "weakest_selected_density": weakest,
            "group_margin": group_margin,
        })
    chromosome_metrics.sort(key=lambda item: (-item["score"], item["chromosome"]))
    best = chromosome_metrics[0]
    second_score = chromosome_metrics[1]["score"] if len(chromosome_metrics) > 1 else 0.0
    chromosome_margin = (
        (best["score"] - second_score) / best["score"] if best["score"] > 0 else 0.0
    )
    total_links = sum(group_links.values())
    selected_fraction = best["selected_links"] / total_links if total_links else 0.0
    return {
        "chromosome": best["chromosome"],
        "groups": best["groups"],
        "chromosome_score": best["score"],
        "second_chromosome_score": second_score,
        "selected_links": best["selected_links"],
        "total_links": total_links,
        "selected_fraction": selected_fraction,
        "chromosome_margin": chromosome_margin,
        "group_margin": best["group_margin"],
        "weakest_selected_density": best["weakest_selected_density"],
    }


def decision_reason(metrics, args):
    if metrics["total_links"] <= 0:
        return "no_assigned_hic"
    if metrics["selected_links"] < args.min_adjusted_links:
        return "insufficient_adjusted_links"
    if metrics["weakest_selected_density"] <= 0:
        return "selected_group_without_hic"
    if metrics["selected_fraction"] < args.min_assigned_fraction:
        return "low_assigned_link_fraction"
    if metrics["chromosome_margin"] < args.min_chromosome_margin:
        return "low_chromosome_margin"
    if metrics["group_margin"] < args.min_group_margin:
        return "low_group_margin"
    return "accepted"


def assign_candidates(
    records, assignments, neighbors, group_re_sites, chromosomes,
    dosage, unsupported, args,
):
    assignments = dict(assignments)
    decisions = {
        unitig: {"round": None, "basis": reason, "metrics": None}
        for unitig, reason in unsupported.items()
    }
    pending = set(records) - set(unsupported)
    round_counts = []
    for round_number in range(1, args.max_rounds + 1):
        accepted = []
        for unitig in sorted(pending, key=lambda item: (-records[item]["length"], item)):
            metrics = score_candidate(
                unitig, assignments, neighbors, group_re_sites,
                chromosomes, args.ploidy, dosage[unitig],
            )
            if decision_reason(metrics, args) == "accepted":
                accepted.append((unitig, metrics))
        if not accepted:
            break
        for unitig, metrics in accepted:
            assignments[unitig] = metrics["groups"]
            decisions[unitig] = {
                "round": round_number,
                "basis": "hic_high_confidence",
                "metrics": metrics,
            }
            for group in metrics["groups"]:
                group_re_sites[group] += records[unitig]["re_sites"] - 1
            pending.remove(unitig)
        round_counts.append(len(accepted))

    if args.low_confidence_policy == "best":
        low_round = args.max_rounds + 1
        while pending:
            accepted = []
            for unitig in sorted(pending, key=lambda item: (-records[item]["length"], item)):
                metrics = score_candidate(
                    unitig, assignments, neighbors, group_re_sites,
                    chromosomes, args.ploidy, dosage[unitig],
                )
                if (
                    metrics["selected_links"] >= args.min_adjusted_links
                    and metrics["weakest_selected_density"] > 0
                ):
                    accepted.append((unitig, metrics))
            if not accepted:
                break
            for unitig, metrics in accepted:
                assignments[unitig] = metrics["groups"]
                decisions[unitig] = {
                    "round": low_round,
                    "basis": "hic_best_available_low_confidence",
                    "metrics": metrics,
                }
                for group in metrics["groups"]:
                    group_re_sites[group] += records[unitig]["re_sites"] - 1
                pending.remove(unitig)
            round_counts.append(len(accepted))
            low_round += 1

    for unitig in sorted(pending):
        metrics = score_candidate(
            unitig, assignments, neighbors, group_re_sites,
            chromosomes, args.ploidy, dosage[unitig],
        )
        decisions[unitig] = {
            "round": None,
            "basis": decision_reason(metrics, args),
            "metrics": metrics,
        }
    return assignments, decisions, round_counts


def copy_validated_fasta(source_path, destination_handle, expected_ids):
    observed = []
    last_line_ended = True
    with source_path.open("rb") as source:
        for line in source:
            destination_handle.write(line)
            last_line_ended = line.endswith(b"\n")
            if line.startswith(b">"):
                observed.append(line[1:].split()[0].decode())
    if not last_line_ended:
        destination_handle.write(b"\n")
    if observed != expected_ids:
        raise ValueError(f"FASTA membership/order mismatch in {source_path}")


def validate_state(state, candidates, candidate_assignments, dosage, assembly_metadata):
    violations = []
    seed_assignments = state["assignments"]
    for unitig, groups in seed_assignments.items():
        if candidate_assignments.get(unitig) != groups:
            violations.append(
                ("seed_changed", unitig, ",".join(groups), str(candidate_assignments.get(unitig)))
            )
    for unitig in candidates:
        groups = candidate_assignments.get(unitig)
        if groups is not None and len(groups) != dosage[unitig]:
            violations.append(("dosage", unitig, str(dosage[unitig]), str(len(groups))))
    for unitig, record in state["inventory"].items():
        if record["status"] == "unassigned" and unitig in candidate_assignments:
            violations.append(("recluster_deferred_changed", unitig, "unassigned", "assigned"))
    inventory_ids = set(state["inventory"]) | set(candidates)
    if inventory_ids != set(assembly_metadata):
        missing = sorted(set(assembly_metadata) - inventory_ids)
        extra = sorted(inventory_ids - set(assembly_metadata))
        for unitig in missing[:20]:
            violations.append(("partition_missing", unitig, "present", "absent"))
        for unitig in extra[:20]:
            violations.append(("partition_extra", unitig, "absent", "present"))
    for unitig, length in assembly_metadata.items():
        if unitig in candidates:
            observed = candidates[unitig]["length"]
        elif unitig in state["inventory"]:
            observed = state["inventory"][unitig]["length"]
        else:
            continue
        if observed != length:
            violations.append(("length", unitig, str(length), str(observed)))
    return violations


def metric_text(metrics, key, digits=6):
    if metrics is None:
        return ""
    return f"{metrics[key]:.{digits}f}"


def write_outputs(
    output_directory, state, candidates, candidate_order, contig_types,
    dosage, dosage_source, assignments, decisions, round_counts,
    assembly_metadata, raw_link_records, relevant_pairs, violations, args,
):
    temporary = Path(tempfile.mkdtemp(prefix=".rescue.", dir=output_directory))
    try:
        group_members = {
            group: list(members) for group, members in state["group_members"].items()
        }
        rescued = [unitig for unitig in candidate_order if unitig in assignments]
        for unitig in rescued:
            for group in assignments[unitig]:
                group_members[group].append(unitig)

        for group in sorted(group_members):
            fasta_path = temporary / f"{group}.reassignment.fa"
            with fasta_path.open("wb") as destination:
                copy_validated_fasta(
                    state["group_fasta_paths"][group], destination,
                    state["group_members"][group],
                )
                for unitig in candidate_order:
                    if group in assignments.get(unitig, ()):
                        destination.write(f">{unitig}\n".encode())
                        sequence = candidates[unitig]["sequence"]
                        for start in range(0, len(sequence), 80):
                            destination.write((sequence[start:start + 80] + "\n").encode())
            with (temporary / f"{group}.reassignment.txt").open("w") as handle:
                handle.writelines(unitig + "\n" for unitig in group_members[group])
            with (temporary / f"{group}.txt").open("w") as handle:
                handle.write("#Contig\tRECounts\tLength\n")
                for unitig in group_members[group]:
                    record = candidates.get(unitig) or state["inventory"][unitig]
                    handle.write(f"{unitig}\t{record['re_sites']}\t{record['length']}\n")

        cluster_text = "".join(
            f"{group}\t{' '.join(group_members[group])}\n"
            for group in sorted(group_members)
        )
        (temporary / "group.reassignment.cluster.txt").write_text(cluster_text)
        (temporary / "merge.group.reassignment.cluster.txt").write_text(cluster_text)

        fieldnames = [
            "unitig", "length", "re_sites", "source", "contig_type", "dosage",
            "dosage_source", "status", "chromosome", "groups", "group_count",
            "assignment_round", "assignment_basis", "decision_adjusted_selected_hic_links",
            "decision_adjusted_total_hic_links", "decision_selected_hic_fraction",
            "decision_chromosome_score", "decision_second_chromosome_score",
            "decision_chromosome_margin", "decision_group_margin",
        ]
        with (temporary / "rescue_assignments.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for unitig in assembly_metadata:
                if unitig in candidates:
                    record = candidates[unitig]
                    decision = decisions[unitig]
                    metrics = decision["metrics"]
                    groups = assignments.get(unitig, ())
                    status = "assigned" if groups else "unassigned"
                    chromosome = (
                        groups[0].rsplit("_group", 1)[0]
                        if groups else (metrics["chromosome"] if metrics else "")
                    )
                    row = {
                        "unitig": unitig,
                        "length": record["length"],
                        "re_sites": record["re_sites"],
                        "source": "rescue_candidate",
                        "contig_type": contig_types.get(unitig, "missing"),
                        "dosage": dosage.get(unitig, ""),
                        "dosage_source": dosage_source[unitig],
                        "status": status,
                        "chromosome": chromosome,
                        "groups": ",".join(groups),
                        "group_count": len(groups),
                        "assignment_round": (
                            decision["round"] if decision["round"] is not None else ""
                        ),
                        "assignment_basis": decision["basis"],
                        "decision_adjusted_selected_hic_links": metric_text(metrics, "selected_links"),
                        "decision_adjusted_total_hic_links": metric_text(metrics, "total_links"),
                        "decision_selected_hic_fraction": metric_text(metrics, "selected_fraction"),
                        "decision_chromosome_score": metric_text(metrics, "chromosome_score", 12),
                        "decision_second_chromosome_score": metric_text(
                            metrics, "second_chromosome_score", 12
                        ),
                        "decision_chromosome_margin": metric_text(metrics, "chromosome_margin"),
                        "decision_group_margin": metric_text(metrics, "group_margin"),
                    }
                else:
                    record = state["inventory"][unitig]
                    groups = assignments.get(unitig, ())
                    row = {
                        "unitig": unitig,
                        "length": record["length"],
                        "re_sites": record["re_sites"],
                        "source": record["source"],
                        "contig_type": record["contig_type"],
                        "dosage": record["dosage"],
                        "dosage_source": record["dosage_source"],
                        "status": record["status"],
                        "chromosome": record["chromosome"],
                        "groups": ",".join(groups),
                        "group_count": len(groups),
                        "assignment_round": "",
                        "assignment_basis": (
                            "fixed_recluster_assignment" if groups
                            else f"deferred_by_recluster:{record['assignment_basis']}"
                        ),
                        "decision_adjusted_selected_hic_links": "",
                        "decision_adjusted_total_hic_links": "",
                        "decision_selected_hic_fraction": "",
                        "decision_chromosome_score": "",
                        "decision_second_chromosome_score": "",
                        "decision_chromosome_margin": "",
                        "decision_group_margin": "",
                    }
                writer.writerow(row)

        final_unassigned = [unitig for unitig in assembly_metadata if unitig not in assignments]
        with (temporary / "unassigned_unitigs.txt").open("w") as handle:
            handle.writelines(unitig + "\n" for unitig in final_unassigned)
        with (temporary / "unassigned_unitigs.fa").open("wb") as destination:
            for chromosome, path in state["unassigned_fasta_paths"]:
                expected = [
                    unitig for unitig, record in state["inventory"].items()
                    if record["chromosome"] == chromosome
                    and record["status"] == "unassigned"
                ]
                copy_validated_fasta(path, destination, expected)
            for unitig in candidate_order:
                if unitig not in assignments:
                    destination.write(f">{unitig}\n".encode())
                    sequence = candidates[unitig]["sequence"]
                    for start in range(0, len(sequence), 80):
                        destination.write((sequence[start:start + 80] + "\n").encode())

        with (temporary / "rescue_validation.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["violation", "unitig", "expected", "observed"])
            writer.writerows(violations)

        assigned_ids = set(assignments)
        final_unassigned = [unitig for unitig in assembly_metadata if unitig not in assignments]
        candidate_unassigned = [unitig for unitig in candidates if unitig not in assignments]
        basis_counts = Counter(decisions[unitig]["basis"] for unitig in candidates)
        group_summary = {
            group: {
                "unitigs": len(members),
                "bp": sum(
                    (candidates.get(unitig) or state["inventory"][unitig])["length"]
                    for unitig in members
                ),
                "re_sites": sum(
                    (candidates.get(unitig) or state["inventory"][unitig])["re_sites"]
                    for unitig in members
                ),
            }
            for group, members in sorted(group_members.items())
        }
        summary = {
            "inputs": {
                "assembly_fasta": str(args.assembly_fasta.resolve()),
                "candidate_fasta": str(args.candidate_fasta.resolve()),
                "contig_type": str(args.contig_type.resolve()),
                "full_links": str(args.full_links.resolve()),
                "recluster_directory": str(args.recluster_dir.resolve()),
            },
            "parameters": {
                "ploidy": args.ploidy,
                "motif": args.RE,
                "flank": args.flank,
                "hic_link_normalization": args.hic_link_normalization,
                "min_adjusted_links": args.min_adjusted_links,
                "min_chromosome_margin": args.min_chromosome_margin,
                "min_group_margin": args.min_group_margin,
                "min_assigned_fraction": args.min_assigned_fraction,
                "max_rounds": args.max_rounds,
                "low_confidence_policy": args.low_confidence_policy,
                "unknown_dosage_policy": args.unknown_dosage_policy,
            },
            "input": {
                "unitigs": len(assembly_metadata),
                "bp": sum(assembly_metadata.values()),
            },
            "recluster_fixed": {
                "unitigs": len(state["assignments"]),
                "bp": sum(
                    state["inventory"][unitig]["length"] for unitig in state["assignments"]
                ),
            },
            "recluster_deferred": {
                "unitigs": sum(
                    record["status"] == "unassigned"
                    for record in state["inventory"].values()
                ),
                "bp": sum(
                    record["length"] for record in state["inventory"].values()
                    if record["status"] == "unassigned"
                ),
            },
            "candidates": {
                "unitigs": len(candidates),
                "bp": sum(record["length"] for record in candidates.values()),
                "rescued_unitigs": len(rescued),
                "rescued_bp": sum(candidates[unitig]["length"] for unitig in rescued),
                "unassigned_unitigs": len(candidate_unassigned),
                "unassigned_bp": sum(
                    candidates[unitig]["length"] for unitig in candidate_unassigned
                ),
            },
            "final": {
                "assigned_unitigs": len(assigned_ids),
                "assigned_bp": sum(assembly_metadata[unitig] for unitig in assigned_ids),
                "unassigned_unitigs": len(final_unassigned),
                "unassigned_bp": sum(
                    assembly_metadata[unitig] for unitig in final_unassigned
                ),
                "memberships": sum(len(groups) for groups in assignments.values()),
            },
            "round_assignment_counts": round_counts,
            "candidate_basis_counts": dict(sorted(basis_counts.items())),
            "candidate_type_counts": dict(sorted(Counter(
                contig_types.get(unitig, "missing") for unitig in candidates
            ).items())),
            "hic": {
                "raw_link_records": raw_link_records,
                "relevant_pairs": relevant_pairs,
            },
            "groups": group_summary,
            "validation": {
                "violations": len(violations),
                "dosage_errors": sum(item[0] == "dosage" for item in violations),
                "changed_recluster_assignments": sum(
                    item[0] == "seed_changed" for item in violations
                ),
                "changed_recluster_deferred": sum(
                    item[0] == "recluster_deferred_changed" for item in violations
                ),
                "partition_errors": sum(
                    item[0].startswith("partition_") for item in violations
                ),
                "length_errors": sum(item[0] == "length" for item in violations),
            },
        }
        with (temporary / "rescue_summary.json").open("w") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")

        generated = {path.name for path in temporary.iterdir()}
        stale_patterns = [
            "*_group*.reassignment.fa", "*_group*.reassignment.txt", "*_group*.txt",
            "group.reassignment.cluster.txt", "merge.group.reassignment.cluster.txt",
            "unassigned_unitigs.txt", "unassigned_unitigs.fa", "rescue_assignments.tsv",
            "rescue_summary.json", "rescue_validation.tsv", "full.links.txt",
            "g1g2g3g4.reassignment.fa", "all_groups.reassignment.fa",
        ]
        for pattern in stale_patterns:
            for old_path in output_directory.glob(pattern):
                if old_path.name not in generated and old_path.is_file():
                    old_path.unlink()
        for path in temporary.iterdir():
            os.replace(path, output_directory / path.name)
        return summary
    finally:
        shutil.rmtree(temporary, ignore_errors=True)


def run(args):
    if args.ploidy < 2:
        raise ValueError("--ploidy must be at least two")
    if args.flank is not None and args.flank < 0:
        raise ValueError("--flank must be non-negative")
    if args.min_adjusted_links < 0:
        raise ValueError("--min-adjusted-links must be non-negative")
    for name in ("min_chromosome_margin", "min_group_margin", "min_assigned_fraction"):
        if not 0 <= getattr(args, name) <= 1:
            raise ValueError(f"--{name.replace('_', '-')} must be between zero and one")
    if args.max_rounds < 1:
        raise ValueError("--max-rounds must be positive")
    if not args.RE:
        raise ValueError("--RE must not be empty")

    output_directory = args.output_dir.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    state = load_recluster_state(args.recluster_dir.resolve(), args.ploidy)
    candidates, candidate_order = parse_fasta(
        args.candidate_fasta, args.RE.upper(), args.flank
    )
    overlap = set(candidates) & set(state["inventory"])
    if overlap:
        raise ValueError(
            f"Rescue candidates already occur in recluster input: {sorted(overlap)[0]}"
        )
    contig_types = parse_contig_types(args.contig_type)
    dosage, dosage_source, unsupported = candidate_dosages(
        candidates, contig_types, args.ploidy, args.unknown_dosage_policy
    )
    all_dosage = {
        unitig: record["dosage"]
        for unitig, record in state["inventory"].items()
        if record["status"] == "assigned"
    }
    all_dosage.update(dosage)
    neighbors, raw_link_records, relevant_pairs = load_relevant_links(
        args.full_links,
        candidates,
        state["assignments"],
        all_dosage,
        args.hic_link_normalization,
    )
    group_re_sites = initial_group_re_sites(
        state["assignments"], state["inventory"], state["group_members"]
    )
    assignments, decisions, round_counts = assign_candidates(
        candidates, state["assignments"], neighbors, group_re_sites,
        state["chromosomes"], dosage, unsupported, args,
    )
    assembly_metadata, assembly_order = scan_fasta_metadata(args.assembly_fasta)
    violations = validate_state(
        state, candidates, assignments, dosage, assembly_metadata
    )
    if violations:
        raise RuntimeError(f"Rescue validation failed with {len(violations)} violations")
    ordered_metadata = {unitig: assembly_metadata[unitig] for unitig in assembly_order}
    summary = write_outputs(
        output_directory, state, candidates, candidate_order, contig_types,
        dosage, dosage_source, assignments, decisions, round_counts,
        ordered_metadata, raw_link_records, relevant_pairs, violations, args,
    )
    print(
        f"Rescued {summary['candidates']['rescued_unitigs']}/"
        f"{summary['candidates']['unitigs']} chromosome-unassigned unitigs; "
        f"final assigned={summary['final']['assigned_unitigs']}/"
        f"{summary['input']['unitigs']}"
    )
    return summary


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Conservative chromosome-and-haplotype rescue for PHap v2"
    )
    parser.add_argument("--assembly-fasta", "--draft_fasta", required=True, type=Path)
    parser.add_argument("--candidate-fasta", "--fasta", required=True, type=Path)
    parser.add_argument("--contig-type", "--contig_type", required=True, type=Path)
    parser.add_argument("--full-links", "--full_links", required=True, type=Path)
    parser.add_argument(
        "--hic-link-normalization",
        choices=["dosage", "raw"],
        default="dosage",
        help=(
            "Normalize each Hi-C count by the product of the two unitig dosages, "
            "or use raw read-pair counts [dosage]"
        ),
    )
    parser.add_argument("--recluster-dir", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--ploidy", type=int, default=4)
    parser.add_argument("--RE", default="GATC")
    parser.add_argument("--flank", type=int)
    parser.add_argument("--min-adjusted-links", type=float, default=5.0)
    parser.add_argument("--min-chromosome-margin", type=float, default=0.10)
    parser.add_argument("--min-group-margin", type=float, default=0.10)
    parser.add_argument("--min-assigned-fraction", type=float, default=0.0)
    parser.add_argument("--max-rounds", type=int, default=10)
    parser.add_argument(
        "--low-confidence-policy", choices=["defer", "best"], default="defer"
    )
    parser.add_argument(
        "--unknown-dosage-policy",
        choices=["defer", "error", "haplotig"], default="defer",
    )
    return parser.parse_args()


def main():
    run(parse_arguments())


if __name__ == "__main__":
    main()
