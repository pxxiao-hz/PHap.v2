#!/usr/bin/env python3
"""Reassign uncertain chromosome unitigs to trusted haplotype groups."""

from __future__ import annotations

import argparse
import csv
import itertools
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


DOSAGE_BY_TYPE = {
    "haplotig": 1,
    "diplotig": 2,
    "triplotig": 3,
    "tetraplotig": 4,
}


def parse_fasta(path: Path, motif: str, flank: Optional[int]):
    records = {}
    order = []
    current = None
    chunks = []

    def finish():
        if current is None:
            return
        sequence = "".join(chunks).upper()
        if flank is None:
            re_sequence = sequence
        else:
            flank_size = min(flank, len(sequence) // 2)
            re_sequence = sequence[:flank_size] + sequence[-flank_size:] if flank_size else ""
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
                    raise ValueError(f"Sequence before first FASTA header at {path}:{line_number}")
                chunks.append(line.strip())
    finish()
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records, order


def parse_contig_types(path: Path, unitigs, unknown_policy: str):
    types = {}
    with path.open() as handle:
        header = handle.readline().split()
        try:
            id_index = header.index("contig_ID")
            type_index = header.index("contig_type")
        except ValueError as exc:
            raise ValueError("contig type file must contain contig_ID and contig_type columns") from exc
        for line_number, line in enumerate(handle, 2):
            fields = line.split()
            if len(fields) <= max(id_index, type_index):
                continue
            unitig = fields[id_index]
            if unitig in types:
                raise ValueError(f"Duplicate contig type for {unitig} at {path}:{line_number}")
            types[unitig] = fields[type_index]

    missing = sorted(set(unitigs) - set(types))
    unsupported = sorted(
        unitig for unitig in unitigs
        if unitig in types and types[unitig] not in DOSAGE_BY_TYPE
    )
    if unknown_policy == "error" and (missing or unsupported):
        examples = missing[:3] + unsupported[:3]
        raise ValueError(
            f"Missing or unsupported contig type for {len(missing) + len(unsupported)} "
            f"chromosome unitigs: {', '.join(examples)}"
        )
    for unitig in missing + unsupported:
        types[unitig] = "haplotig"
    inferred = set(missing + unsupported)
    return {unitig: types[unitig] for unitig in unitigs}, inferred


def parse_seed_clusters(path: Path, records, dosage, ploidy: int):
    assignments = {}
    seen_groups = set()
    chromosome = None
    pattern = re.compile(r"^(?:.*_)?group([0-9]+)$")
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            match = pattern.match(fields[0])
            if not match:
                raise ValueError(f"Invalid seed group name {fields[0]} at {path}:{line_number}")
            group = int(match.group(1)) - 1
            if group < 0 or group >= ploidy or group in seen_groups:
                raise ValueError(f"Duplicate or out-of-range group at {path}:{line_number}")
            seen_groups.add(group)
            if "_group" in fields[0]:
                chromosome = fields[0].rsplit("_group", 1)[0]
                members = fields[1:]
            else:
                if len(fields) < 2:
                    raise ValueError(f"Malformed seed group at {path}:{line_number}")
                try:
                    declared = int(fields[1])
                except ValueError as exc:
                    raise ValueError(f"Missing seed member count at {path}:{line_number}") from exc
                members = fields[2:]
                if declared != len(members):
                    raise ValueError(
                        f"Seed member count mismatch at {path}:{line_number}: "
                        f"declared {declared}, observed {len(members)}"
                    )
            if len(members) != len(set(members)):
                raise ValueError(f"Duplicate unitig in group at {path}:{line_number}")
            for unitig in members:
                if unitig not in records:
                    raise ValueError(f"Seed unitig {unitig} is absent from chromosome FASTA")
                assignments.setdefault(unitig, set()).add(group)

    expected_groups = set(range(ploidy))
    if seen_groups != expected_groups:
        missing = sorted(group + 1 for group in expected_groups - seen_groups)
        raise ValueError(f"Seed cluster file is missing groups: {missing}")
    for unitig, groups in assignments.items():
        if len(groups) != dosage[unitig]:
            raise ValueError(
                f"Seed dosage mismatch for {unitig}: expected {dosage[unitig]}, observed {len(groups)}"
            )
    if not assignments:
        raise ValueError(f"No seed unitigs found in {path}")
    return {unitig: tuple(sorted(groups)) for unitig, groups in assignments.items()}, chromosome


def parse_cluster_evidence(
    path, seed_assignments, trusted_bases, min_links, min_margin, min_fraction
):
    evidence = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "unitig", "groups", "assignment_basis",
            "adjusted_assigned_hic_links", "assigned_hic_fraction",
            "hic_density_margin",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"Missing columns in {path}: {', '.join(sorted(missing))}"
            )
        for line_number, row in enumerate(reader, 2):
            unitig = row["unitig"]
            if unitig in evidence:
                raise ValueError(f"Duplicate cluster evidence for {unitig} at {path}:{line_number}")
            if unitig not in seed_assignments:
                raise ValueError(
                    f"Cluster evidence unitig {unitig} is absent from the seed groups"
                )
            try:
                groups = tuple(sorted(int(group) - 1 for group in row["groups"].split(",") if group))
            except ValueError as exc:
                raise ValueError(f"Invalid cluster groups for {unitig} at {path}:{line_number}") from exc
            if groups != seed_assignments[unitig]:
                raise ValueError(
                    f"Cluster evidence/group mismatch for {unitig}: "
                    f"evidence={groups}, groups={seed_assignments[unitig]}"
                )
            basis = row["assignment_basis"]
            try:
                assigned_links = float(row["adjusted_assigned_hic_links"])
                assigned_fraction = float(row["assigned_hic_fraction"])
                margin = float(row["hic_density_margin"])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid cluster evidence metrics for {unitig} at {path}:{line_number}"
                ) from exc
            trusted = (
                basis in trusted_bases
                and assigned_links >= min_links
                and assigned_fraction >= min_fraction
                and margin >= min_margin
            )
            evidence[unitig] = {
                "assignment_basis": basis,
                "assigned_links": assigned_links,
                "assigned_fraction": assigned_fraction,
                "margin": margin,
                "trusted": trusted,
            }
    missing_seeds = sorted(set(seed_assignments) - set(evidence))
    if missing_seeds:
        raise ValueError(
            f"Cluster evidence is missing {len(missing_seeds)} seed unitigs: {missing_seeds[0]}"
        )
    return evidence


def load_relevant_links(path: Path, units, dosage, normalization="dosage"):
    if normalization not in {"dosage", "raw"}:
        raise ValueError(
            f"Unsupported Hi-C link normalization: {normalization}; expected dosage or raw"
        )
    with path.open("rb") as handle:
        raw_links = pickle.load(handle)
    if not isinstance(raw_links, dict):
        raise ValueError("Hi-C link pickle must contain a dictionary")
    links = defaultdict(float)
    invalid = 0
    unit_set = set(units)
    for key, value in raw_links.items():
        if not isinstance(key, tuple) or len(key) != 2:
            invalid += 1
            continue
        unitig1, unitig2 = key
        if unitig1 == unitig2 or unitig1 not in unit_set or unitig2 not in unit_set:
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
        links[pair] += count
    if invalid:
        raise ValueError(f"Hi-C link pickle contains {invalid} invalid records")
    neighbors = defaultdict(dict)
    for (unitig1, unitig2), count in links.items():
        neighbors[unitig1][unitig2] = count
        neighbors[unitig2][unitig1] = count
    return dict(neighbors), len(raw_links), len(links)


def group_re_sites(assignments, records, ploidy, exclude=None):
    totals = [1] * ploidy
    for unitig, groups in assignments.items():
        if unitig == exclude:
            continue
        for group in groups:
            totals[group] += records[unitig]["re_sites"] - 1
    return totals


def score_unitig(
    unitig, assignments, neighbors, records, dosage, ploidy, forced_groups=None
):
    group_links = [0.0] * ploidy
    for neighbor, count in neighbors.get(unitig, {}).items():
        for group in assignments.get(neighbor, ()):
            group_links[group] += count
    re_totals = group_re_sites(assignments, records, ploidy, exclude=unitig)
    densities = [group_links[group] / max(re_totals[group], 1) for group in range(ploidy)]
    ranked = sorted(range(ploidy), key=lambda group: (-densities[group], group))
    selected = (
        tuple(sorted(forced_groups))
        if forced_groups is not None
        else tuple(sorted(ranked[: dosage[unitig]]))
    )
    selected_links = sum(group_links[group] for group in selected)
    total_links = sum(group_links)
    selected_fraction = selected_links / total_links if total_links else 0.0
    weakest_group = min(selected, key=lambda group: (densities[group], group), default=None)
    weakest_selected = densities[weakest_group] if weakest_group is not None else 0.0
    unselected = [group for group in range(ploidy) if group not in selected]
    strongest_group = max(
        unselected, key=lambda group: (densities[group], -group), default=None
    )
    strongest_unselected = (
        densities[strongest_group] if strongest_group is not None else 0.0
    )
    if not unselected:
        margin = 1.0
    elif weakest_selected > 0:
        margin = (weakest_selected - strongest_unselected) / weakest_selected
    else:
        margin = 0.0
    return {
        "groups": selected,
        "group_links": group_links,
        "densities": densities,
        "group_re_sites": re_totals,
        "selected_links": selected_links,
        "total_links": total_links,
        "selected_fraction": selected_fraction,
        "margin": margin,
        "weakest_selected_density": weakest_selected,
        "weakest_selected_group": weakest_group,
        "strongest_unselected_group": strongest_group,
    }


def decision_reason(
    metrics, min_links, min_margin, min_fraction
):
    if metrics["total_links"] <= 0:
        return "no_assigned_hic"
    if metrics["selected_links"] < min_links:
        return "insufficient_adjusted_links"
    if metrics["weakest_selected_density"] <= 0:
        return "selected_group_without_hic"
    if metrics["selected_fraction"] < min_fraction:
        return "low_assigned_link_fraction"
    if metrics["margin"] < min_margin:
        return "low_group_margin"
    return "accepted"


def assign_nonseed_unitigs(
    records, fixed_seed_assignments, original_seed_assignments, neighbors, dosage,
    ploidy, min_links, min_margin, min_fraction, max_rounds, low_confidence_policy,
    reviewed_seed_fallback,
):
    assignments = dict(fixed_seed_assignments)
    reviewed_seeds = set(original_seed_assignments) - set(fixed_seed_assignments)
    decisions = {
        unitig: {"round": 0, "basis": "fixed_seed", "metrics": None}
        for unitig in fixed_seed_assignments
    }
    pending = set(records) - set(assignments)
    all_groups = tuple(range(ploidy))
    for unitig in sorted(tuple(pending)):
        if dosage[unitig] == ploidy:
            assignments[unitig] = all_groups
            basis = (
                "reviewed_seed_dosage_all_groups"
                if unitig in reviewed_seeds else "dosage_all_groups"
            )
            decisions[unitig] = {"round": 0, "basis": basis, "metrics": None}
            pending.remove(unitig)

    round_counts = []
    for round_number in range(1, max_rounds + 1):
        accepted = []
        for unitig in sorted(pending, key=lambda item: (-records[item]["length"], item)):
            metrics = score_unitig(unitig, assignments, neighbors, records, dosage, ploidy)
            if decision_reason(metrics, min_links, min_margin, min_fraction) == "accepted":
                accepted.append((unitig, metrics))
        if not accepted:
            break
        for unitig, metrics in accepted:
            assignments[unitig] = metrics["groups"]
            decisions[unitig] = {
                "round": round_number,
                "basis": (
                    "reviewed_seed_hic_high_confidence"
                    if unitig in reviewed_seeds else "hic_high_confidence"
                ),
                "metrics": metrics,
            }
            pending.remove(unitig)
        round_counts.append(len(accepted))

    if low_confidence_policy == "best":
        low_round = max_rounds + 1
        while True:
            accepted = []
            for unitig in sorted(pending, key=lambda item: (-records[item]["length"], item)):
                metrics = score_unitig(unitig, assignments, neighbors, records, dosage, ploidy)
                if metrics["selected_links"] >= min_links and metrics["weakest_selected_density"] > 0:
                    accepted.append((unitig, metrics))
            if not accepted:
                break
            for unitig, metrics in accepted:
                assignments[unitig] = metrics["groups"]
                decisions[unitig] = {
                    "round": low_round,
                    "basis": (
                        "reviewed_seed_hic_best_available_low_confidence"
                        if unitig in reviewed_seeds else "hic_best_available_low_confidence"
                    ),
                    "metrics": metrics,
                }
                pending.remove(unitig)
            round_counts.append(len(accepted))
            low_round += 1

    for unitig in sorted(pending):
        metrics = score_unitig(unitig, assignments, neighbors, records, dosage, ploidy)
        basis = decision_reason(metrics, min_links, min_margin, min_fraction)
        if unitig in reviewed_seeds:
            if reviewed_seed_fallback == "retain":
                assignments[unitig] = original_seed_assignments[unitig]
                basis = f"reviewed_seed_retained_{basis}"
            else:
                basis = f"reviewed_seed_{basis}"
        decisions[unitig] = {
            "round": None,
            "basis": basis,
            "metrics": metrics,
        }
    return assignments, decisions, round_counts


def assignment_signature(assignments):
    return tuple(sorted((unitig, tuple(groups)) for unitig, groups in assignments.items()))


def parse_relaxed_constraints(path: Path):
    relaxed = set()
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"unitig1", "unitig2"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError(f"Relaxed-constraint table lacks unitig1/unitig2: {path}")
        for row in reader:
            relaxed.add(tuple(sorted((row["unitig1"], row["unitig2"]))))
    return relaxed


def parse_allelic_blocks(path: Path, records, relaxed_pairs=None):
    """Read chromosome-local allelic rows as small, auditable joint blocks."""
    blocks = []
    relaxed_pairs = relaxed_pairs or set()
    seen = set()
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid allelic-table coordinates at {path}:{line_number}"
                ) from exc
            units = tuple(dict.fromkeys(unitig for unitig in fields[3:] if unitig in records))
            if len(units) < 2:
                continue
            if any(
                tuple(sorted(pair)) in relaxed_pairs
                for pair in itertools.combinations(units, 2)
            ):
                continue
            key = frozenset(units)
            if key in seen:
                continue
            seen.add(key)
            blocks.append(
                {
                    "block_id": f"row{line_number}",
                    "chromosome": fields[0],
                    "start": start,
                    "end": end,
                    "unitigs": units,
                }
            )
    return blocks


def enumerate_block_configurations(
    movable, occupied_groups, dosage, ploidy, max_configurations
):
    domains = {
        unitig: tuple(itertools.combinations(range(ploidy), dosage[unitig]))
        for unitig in movable
    }
    configurations = []

    def search(index, used, current):
        if len(configurations) > max_configurations:
            return
        if index == len(movable):
            configurations.append(dict(current))
            return
        unitig = movable[index]
        for groups in domains[unitig]:
            group_set = set(groups)
            if group_set & used:
                continue
            current[unitig] = tuple(groups)
            search(index + 1, used | group_set, current)
            current.pop(unitig, None)

    search(0, set(occupied_groups), {})
    if len(configurations) > max_configurations:
        return []
    return configurations


def joint_reassign_allelic_blocks(
    records,
    assignments,
    decisions,
    fixed_seed_assignments,
    neighbors,
    dosage,
    blocks,
    args,
):
    """Resolve weak unitigs jointly inside strict allelic-table rows.

    High-confidence assignments act as anchors. Low-confidence or retained
    seeds are optimized together under exact dosage and mutual exclusion. A
    block is accepted only by a unique constraint-forced completion or by the
    configured block-margin threshold.
    """
    assignments = dict(assignments)
    locked = set()
    events = []
    low_confidence_prefixes = (
        "reviewed_seed_retained_",
        "reviewed_seed_low_",
        "low_",
        "insufficient_",
        "no_assigned_",
        "selected_group_",
    )

    def is_movable(unitig):
        if unitig in fixed_seed_assignments or unitig in locked:
            return False
        if unitig not in assignments:
            return True
        basis = decisions[unitig]["basis"]
        return basis.startswith(low_confidence_prefixes)

    ordered_blocks = sorted(
        blocks,
        key=lambda block: (
            -sum(unitig in assignments for unitig in block["unitigs"]),
            block["start"],
            block["end"],
            block["block_id"],
        ),
    )
    for block in ordered_blocks:
        units = tuple(block["unitigs"])
        if sum(dosage[unitig] for unitig in units) > args.ploidy:
            continue
        movable = tuple(
            sorted(
                (unitig for unitig in units if is_movable(unitig)),
                key=lambda unitig: (-records[unitig]["length"], unitig),
            )
        )
        if not movable:
            continue
        anchors = tuple(
            unitig for unitig in units
            if unitig not in movable and unitig in assignments
        )
        occupied = set()
        incompatible_anchors = False
        for unitig in anchors:
            groups = set(assignments[unitig])
            if groups & occupied:
                incompatible_anchors = True
                break
            occupied.update(groups)
        if incompatible_anchors or len(anchors) < args.min_allelic_block_anchors:
            continue

        configurations = enumerate_block_configurations(
            movable,
            occupied,
            dosage,
            args.ploidy,
            args.max_allelic_block_configurations,
        )
        if not configurations:
            continue
        support_assignments = {
            unitig: groups for unitig, groups in assignments.items()
            if unitig not in movable
        }
        scored = []
        for configuration in configurations:
            unit_metrics = {
                unitig: score_unitig(
                    unitig,
                    support_assignments,
                    neighbors,
                    records,
                    dosage,
                    args.ploidy,
                    forced_groups=groups,
                )
                for unitig, groups in configuration.items()
            }
            score = sum(
                sum(metrics["densities"][group] for group in metrics["groups"])
                for metrics in unit_metrics.values()
            )
            selected_links = sum(
                metrics["selected_links"] for metrics in unit_metrics.values()
            )
            signature = tuple(
                (unitig, configuration[unitig]) for unitig in movable
            )
            scored.append((score, selected_links, signature, configuration, unit_metrics))
        scored.sort(key=lambda item: (-item[0], -item[1], item[2]))
        best = scored[0]
        second = scored[1] if len(scored) > 1 else None
        acceptance_basis = "none"
        block_margin = 1.0 if second is None else (
            (best[0] - second[0]) / best[0] if best[0] > 0 else 0.0
        )
        if second is None:
            accepted = True
            acceptance_basis = "constraint_forced"
        else:
            accepted = (
                block_margin >= args.min_group_margin
                and best[1] >= args.min_adjusted_links
            )
            if accepted:
                acceptance_basis = "margin"
        if not accepted:
            continue

        changed = []
        for unitig in movable:
            old_groups = assignments.get(unitig, ())
            new_groups = best[3][unitig]
            if old_groups == new_groups:
                continue
            assignments[unitig] = new_groups
            metrics = best[4][unitig]
            metrics.update(
                {
                    "allelic_block_id": block["block_id"],
                    "allelic_block_margin": block_margin,
                    "allelic_block_acceptance_basis": acceptance_basis,
                }
            )
            reviewed = unitig in decisions and decisions[unitig]["basis"].startswith(
                "reviewed_seed_"
            )
            decisions[unitig] = {
                "round": decisions.get(unitig, {}).get("round"),
                "basis": (
                    "reviewed_seed_allelic_block_joint"
                    if reviewed else "allelic_block_joint"
                ),
                "metrics": metrics,
            }
            changed.append((unitig, old_groups, new_groups))
        if not changed:
            continue
        locked.update(movable)
        events.append(
            {
                **block,
                "anchors": anchors,
                "movable": movable,
                "changes": changed,
                "configurations": len(scored),
                "best_score": best[0],
                "second_score": second[0] if second is not None else None,
                "block_margin": block_margin,
                "acceptance_basis": acceptance_basis,
            }
        )
    return assignments, decisions, events


def refine_assignments(
    records, assignments, decisions, fixed_seed_assignments, neighbors, dosage,
    ploidy, min_links, min_margin, min_fraction, max_rounds,
):
    assignments = dict(assignments)
    excluded_support = {
        unitig for unitig, decision in decisions.items()
        if decision["basis"].startswith("reviewed_seed_retained_")
    }
    seen = {assignment_signature(assignments)}
    change_counts = []
    events = []
    stable = False
    oscillation = False

    for refinement_round in range(1, max_rounds + 1):
        support_assignments = {
            unitig: groups for unitig, groups in assignments.items()
            if unitig not in excluded_support
        }
        proposals = []
        for unitig in sorted(
            set(records) - set(fixed_seed_assignments),
            key=lambda item: (-records[item]["length"], item),
        ):
            metrics = score_unitig(
                unitig, support_assignments, neighbors, records, dosage, ploidy
            )
            if decision_reason(metrics, min_links, min_margin, min_fraction) != "accepted":
                continue
            old_groups = assignments.get(unitig, ())
            if metrics["groups"] != old_groups:
                proposals.append((unitig, old_groups, metrics))

        if not proposals:
            change_counts.append(0)
            stable = True
            break

        candidate_assignments = dict(assignments)
        for unitig, _, metrics in proposals:
            candidate_assignments[unitig] = metrics["groups"]
        signature = assignment_signature(candidate_assignments)
        if signature in seen:
            oscillation = True
            break

        assignments = candidate_assignments
        seen.add(signature)
        change_counts.append(len(proposals))
        for unitig, old_groups, metrics in proposals:
            excluded_support.discard(unitig)
            reviewed_seed = decisions[unitig]["basis"].startswith("reviewed_seed_")
            decisions[unitig] = {
                "round": decisions[unitig].get("round"),
                "refinement_round": refinement_round,
                "basis": (
                    "reviewed_seed_hic_refined"
                    if reviewed_seed else "hic_refined"
                ),
                "metrics": metrics,
            }
            events.append(
                {
                    "refinement_round": refinement_round,
                    "unitig": unitig,
                    "length": records[unitig]["length"],
                    "old_groups": old_groups,
                    "new_groups": metrics["groups"],
                    "selected_links": metrics["selected_links"],
                    "total_links": metrics["total_links"],
                    "selected_fraction": metrics["selected_fraction"],
                    "group_margin": metrics["margin"],
                }
            )

    return assignments, decisions, change_counts, events, stable, oscillation


def validate(records, fixed_seed_assignments, assignments, dosage, block_events=()):
    violations = []
    for unitig, groups in assignments.items():
        if unitig not in records:
            violations.append(("assigned_not_in_fasta", unitig, "present", "absent"))
            continue
        if len(groups) != dosage[unitig]:
            violations.append(("dosage", unitig, str(dosage[unitig]), str(len(groups))))
        if len(groups) != len(set(groups)):
            violations.append(("duplicate_group", unitig, "unique", str(groups)))
    for unitig, groups in fixed_seed_assignments.items():
        if assignments.get(unitig) != groups:
            violations.append(("seed_changed", unitig, str(groups), str(assignments.get(unitig))))
    for event in block_events:
        units = event["unitigs"]
        for unitig1, unitig2 in itertools.combinations(units, 2):
            shared = set(assignments.get(unitig1, ())) & set(assignments.get(unitig2, ()))
            if shared:
                violations.append(
                    (
                        "allelic_block_conflict",
                        f"{unitig1},{unitig2}",
                        "disjoint",
                        ",".join(str(group + 1) for group in sorted(shared)),
                    )
                )
    return violations


def write_fasta(handle, unitig, sequence, prefix=None):
    identifier = f"{prefix}_{unitig}" if prefix else unitig
    handle.write(f">{identifier}\n{sequence}\n")


def write_outputs(
    output_directory, chromosome, records, order, contig_types, dosage,
    original_seed_assignments, fixed_seed_assignments, cluster_evidence,
    assignments, decisions, round_counts, refinement_change_counts,
    refinement_events, refinement_stable, refinement_oscillation, neighbors,
    block_events, raw_link_records, relevant_pairs, violations, args,
):
    temporary = Path(tempfile.mkdtemp(prefix=".recluster.", dir=output_directory))
    try:
        group_members = [
            sorted(unitig for unitig, groups in assignments.items() if group in groups)
            for group in range(args.ploidy)
        ]
        for group, members in enumerate(group_members, 1):
            with (temporary / f"group{group}.reassignment.txt").open("w") as handle:
                handle.writelines(unitig + "\n" for unitig in members)
            with (temporary / f"group{group}.reassignment.fa").open("w") as handle:
                for unitig in members:
                    write_fasta(handle, unitig, records[unitig]["sequence"])
            with (temporary / f"group{group}.txt").open("w") as handle:
                handle.write("#Contig\tRECounts\tLength\n")
                for unitig in members:
                    handle.write(f"{unitig}\t{records[unitig]['re_sites']}\t{records[unitig]['length']}\n")

        with (temporary / "group.reassignment.cluster.txt").open("w") as handle:
            for group, members in enumerate(group_members, 1):
                handle.write(f"{chromosome}_group{group}\t{' '.join(members)}\n")

        combined_names = (
            ["g1g2g3g4.reassignment.fa"]
            if args.ploidy == 4
            else ["all_groups.reassignment.fa"]
        )
        for name in combined_names:
            with (temporary / name).open("w") as handle:
                for group, members in enumerate(group_members, 1):
                    for unitig in members:
                        write_fasta(handle, unitig, records[unitig]["sequence"], prefix=f"g{group}")

        fields = [
            "unitig", "length", "re_sites", "contig_type", "dosage", "dosage_source", "status",
            "groups", "group_count", "seed", "original_seed_groups", "cluster_assignment_basis",
            "seed_review_status", "seed_changed", "assignment_round", "refinement_round",
            "assignment_basis",
            "decision_adjusted_assigned_hic_links", "decision_adjusted_total_hic_links",
            "decision_assigned_hic_fraction", "decision_group_margin",
            "allelic_block_id", "allelic_block_margin",
            "allelic_block_acceptance_basis",
            "final_adjusted_assigned_hic_links", "final_adjusted_total_hic_links",
            "final_assigned_hic_fraction", "final_group_margin",
            *[f"g{group + 1}_links" for group in range(args.ploidy)],
            *[f"g{group + 1}_density" for group in range(args.ploidy)],
        ]
        basis_counts = Counter()
        with (temporary / "recluster_assignments.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for unitig in order:
                groups = assignments.get(unitig, ())
                decision_metrics = decisions[unitig]["metrics"]
                final_metrics = score_unitig(
                    unitig,
                    assignments,
                    neighbors,
                    records,
                    dosage,
                    args.ploidy,
                    forced_groups=groups if groups else None,
                )
                basis = decisions[unitig]["basis"]
                basis_counts[basis] += 1
                original_groups = original_seed_assignments.get(unitig)
                if original_groups is None:
                    review_status = "not_seed"
                    changed = ""
                elif unitig in fixed_seed_assignments:
                    review_status = "fixed_trusted"
                    changed = "no"
                elif not groups:
                    review_status = "reviewed_deferred"
                    changed = ""
                elif tuple(groups) == original_groups:
                    review_status = "reviewed_retained"
                    changed = "no"
                else:
                    review_status = "reviewed_changed"
                    changed = "yes"
                row = {
                    "unitig": unitig,
                    "length": records[unitig]["length"],
                    "re_sites": records[unitig]["re_sites"],
                    "contig_type": contig_types[unitig],
                    "dosage": dosage[unitig],
                    "dosage_source": (
                        "default_haplotig" if unitig in args.inferred_dosage else "provided"
                    ),
                    "status": "assigned" if groups else "unassigned",
                    "groups": ",".join(str(group + 1) for group in groups),
                    "group_count": len(groups),
                    "seed": "yes" if unitig in original_seed_assignments else "no",
                    "original_seed_groups": (
                        ",".join(str(group + 1) for group in original_groups)
                        if original_groups is not None else ""
                    ),
                    "cluster_assignment_basis": (
                        cluster_evidence[unitig]["assignment_basis"]
                        if unitig in cluster_evidence else ""
                    ),
                    "seed_review_status": review_status,
                    "seed_changed": changed,
                    "assignment_round": decisions[unitig]["round"] if decisions[unitig]["round"] is not None else "",
                    "refinement_round": decisions[unitig].get("refinement_round", ""),
                    "assignment_basis": basis,
                    "decision_adjusted_assigned_hic_links": (
                        f'{decision_metrics["selected_links"]:.6f}' if decision_metrics else ""
                    ),
                    "decision_adjusted_total_hic_links": (
                        f'{decision_metrics["total_links"]:.6f}' if decision_metrics else ""
                    ),
                    "decision_assigned_hic_fraction": (
                        f'{decision_metrics["selected_fraction"]:.6f}' if decision_metrics else ""
                    ),
                    "decision_group_margin": (
                        f'{decision_metrics["margin"]:.6f}' if decision_metrics else ""
                    ),
                    "allelic_block_id": (
                        decision_metrics.get("allelic_block_id", "")
                        if decision_metrics else ""
                    ),
                    "allelic_block_margin": (
                        f'{decision_metrics["allelic_block_margin"]:.6f}'
                        if decision_metrics and "allelic_block_margin" in decision_metrics else ""
                    ),
                    "allelic_block_acceptance_basis": (
                        decision_metrics.get("allelic_block_acceptance_basis", "")
                        if decision_metrics else ""
                    ),
                    "final_adjusted_assigned_hic_links": f'{final_metrics["selected_links"]:.6f}',
                    "final_adjusted_total_hic_links": f'{final_metrics["total_links"]:.6f}',
                    "final_assigned_hic_fraction": f'{final_metrics["selected_fraction"]:.6f}',
                    "final_group_margin": f'{final_metrics["margin"]:.6f}',
                }
                row.update({f"g{group + 1}_links": f'{final_metrics["group_links"][group]:.6f}' for group in range(args.ploidy)})
                row.update({f"g{group + 1}_density": f'{final_metrics["densities"][group]:.12f}' for group in range(args.ploidy)})
                writer.writerow(row)

        with (temporary / "recluster_refinement.tsv").open("w", newline="") as handle:
            fields = [
                "refinement_round", "unitig", "length", "old_groups", "new_groups",
                "adjusted_selected_hic_links", "adjusted_total_hic_links",
                "selected_hic_fraction", "group_margin",
            ]
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for event in refinement_events:
                writer.writerow(
                    {
                        "refinement_round": event["refinement_round"],
                        "unitig": event["unitig"],
                        "length": event["length"],
                        "old_groups": ",".join(
                            str(group + 1) for group in event["old_groups"]
                        ),
                        "new_groups": ",".join(
                            str(group + 1) for group in event["new_groups"]
                        ),
                        "adjusted_selected_hic_links": f'{event["selected_links"]:.6f}',
                        "adjusted_total_hic_links": f'{event["total_links"]:.6f}',
                        "selected_hic_fraction": f'{event["selected_fraction"]:.6f}',
                        "group_margin": f'{event["group_margin"]:.6f}',
                    }
                )

        with (temporary / "allelic_block_reassignments.tsv").open(
            "w", newline=""
        ) as handle:
            fields = [
                "block_id", "chromosome", "start", "end", "unitigs", "anchors",
                "movable", "changed_unitigs", "old_groups", "new_groups",
                "configurations", "best_score", "second_score", "block_margin",
                "acceptance_basis",
            ]
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for event in block_events:
                writer.writerow(
                    {
                        "block_id": event["block_id"],
                        "chromosome": event["chromosome"],
                        "start": event["start"],
                        "end": event["end"],
                        "unitigs": ",".join(event["unitigs"]),
                        "anchors": ",".join(event["anchors"]),
                        "movable": ",".join(event["movable"]),
                        "changed_unitigs": ",".join(
                            change[0] for change in event["changes"]
                        ),
                        "old_groups": ";".join(
                            f'{unitig}:{",".join(str(group + 1) for group in old)}'
                            for unitig, old, _ in event["changes"]
                        ),
                        "new_groups": ";".join(
                            f'{unitig}:{",".join(str(group + 1) for group in new)}'
                            for unitig, _, new in event["changes"]
                        ),
                        "configurations": event["configurations"],
                        "best_score": f'{event["best_score"]:.12g}',
                        "second_score": (
                            f'{event["second_score"]:.12g}'
                            if event["second_score"] is not None else ""
                        ),
                        "block_margin": f'{event["block_margin"]:.6f}',
                        "acceptance_basis": event["acceptance_basis"],
                    }
                )

        unassigned = [unitig for unitig in order if unitig not in assignments]
        with (temporary / "unassigned_unitigs.txt").open("w") as handle:
            handle.writelines(unitig + "\n" for unitig in unassigned)
        with (temporary / "unassigned_unitigs.fa").open("w") as handle:
            for unitig in unassigned:
                write_fasta(handle, unitig, records[unitig]["sequence"])
        with (temporary / "recluster_validation.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["violation", "unitig", "expected", "observed"])
            writer.writerows(violations)

        assigned_bp = sum(records[unitig]["length"] for unitig in assignments)
        input_bp = sum(item["length"] for item in records.values())
        seed_bp = sum(records[unitig]["length"] for unitig in original_seed_assignments)
        fixed_seed_bp = sum(records[unitig]["length"] for unitig in fixed_seed_assignments)
        reviewed_seeds = set(original_seed_assignments) - set(fixed_seed_assignments)
        reviewed_assigned = reviewed_seeds & set(assignments)
        reviewed_changed = {
            unitig for unitig in reviewed_assigned
            if assignments[unitig] != original_seed_assignments[unitig]
        }
        nonseed_assigned = set(assignments) - set(original_seed_assignments)
        groups_summary = {
            f"g{group}": {
                "unitigs": len(members),
                "bp": sum(records[unitig]["length"] for unitig in members),
                "re_sites": sum(records[unitig]["re_sites"] for unitig in members),
            }
            for group, members in enumerate(group_members, 1)
        }
        summary = {
            "chromosome": chromosome,
            "inputs": {
                "fasta": str(args.fasta.resolve()),
                "contig_type": str(args.contig_type.resolve()),
                "full_links": str(args.full_links.resolve()),
                "clusters_file": str(args.clusters_file.resolve()),
                "cluster_assignments": (
                    str(args.cluster_assignments.resolve())
                    if args.cluster_assignments else None
                ),
                "allelic_table": (
                    str(args.allelic_table.resolve()) if args.allelic_table else None
                ),
                "relaxed_constraints": (
                    str(args.relaxed_constraints.resolve())
                    if args.relaxed_constraints else None
                ),
            },
            "parameters": {
                "ploidy": args.ploidy,
                "motif": args.RE,
                "flank": args.flank,
                "hic_link_normalization": args.hic_link_normalization,
                "min_adjusted_links": args.min_adjusted_links,
                "min_group_margin": args.min_group_margin,
                "min_allelic_block_anchors": args.min_allelic_block_anchors,
                "max_allelic_block_configurations": args.max_allelic_block_configurations,
                "min_assigned_fraction": args.min_assigned_fraction,
                "max_rounds": args.max_rounds,
                "refinement_rounds": args.refinement_rounds,
                "low_confidence_policy": args.low_confidence_policy,
                "unknown_dosage_policy": args.unknown_dosage_policy,
                "reviewed_seed_fallback": args.reviewed_seed_fallback,
                "trusted_seed_bases": sorted(args.trusted_seed_bases),
            },
            "inferred_haplotig_unitigs": len(args.inferred_dosage),
            "input": {"unitigs": len(records), "bp": input_bp},
            "seed": {"unitigs": len(original_seed_assignments), "bp": seed_bp},
            "fixed_seed": {"unitigs": len(fixed_seed_assignments), "bp": fixed_seed_bp},
            "reviewed_seed": {
                "unitigs": len(reviewed_seeds),
                "bp": sum(records[unitig]["length"] for unitig in reviewed_seeds),
                "assigned_unitigs": len(reviewed_assigned),
                "assigned_bp": sum(records[unitig]["length"] for unitig in reviewed_assigned),
                "changed_unitigs": len(reviewed_changed),
                "changed_bp": sum(records[unitig]["length"] for unitig in reviewed_changed),
                "deferred_unitigs": len(reviewed_seeds - set(assignments)),
                "deferred_bp": sum(
                    records[unitig]["length"] for unitig in reviewed_seeds - set(assignments)
                ),
            },
            "assigned": {
                "unitigs": len(assignments), "bp": assigned_bp,
                "new_unitigs": len(nonseed_assigned),
                "new_bp": sum(records[unitig]["length"] for unitig in nonseed_assigned),
            },
            "unassigned": {"unitigs": len(unassigned), "bp": input_bp - assigned_bp},
            "round_assignment_counts": round_counts,
            "refinement_round_change_counts": refinement_change_counts,
            "refinement": {
                "events": len(refinement_events),
                "stable": refinement_stable,
                "oscillation": refinement_oscillation,
            },
            "allelic_block_joint": {
                "accepted_blocks": len(block_events),
                "changed_unitigs": sum(
                    len(event["changes"]) for event in block_events
                ),
                "acceptance_basis_counts": dict(
                    sorted(Counter(event["acceptance_basis"] for event in block_events).items())
                ),
            },
            "assignment_basis_counts": dict(sorted(basis_counts.items())),
            "hic": {"raw_link_records": raw_link_records, "chromosome_relevant_pairs": relevant_pairs},
            "groups": groups_summary,
            "validation": {
                "violations": len(violations),
                "dosage_errors": sum(item[0] == "dosage" for item in violations),
                "changed_fixed_seeds": sum(item[0] == "seed_changed" for item in violations),
                "allelic_block_conflicts": sum(
                    item[0] == "allelic_block_conflict" for item in violations
                ),
                "partition_errors": len(set(records) - set(assignments) - set(unassigned)),
            },
        }
        with (temporary / "recluster_summary.json").open("w") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")

        obsolete_names = ["full.links.txt", "flank.links.txt"]
        obsolete_names.append(
            "all_groups.reassignment.fa"
            if args.ploidy == 4
            else "g1g2g3g4.reassignment.fa"
        )
        for name in obsolete_names:
            old_path = output_directory / name
            if old_path.exists():
                old_path.unlink()
        generated = {path.name for path in temporary.iterdir()}
        for pattern in ("group*.reassignment.fa", "group*.reassignment.txt", "group[0-9]*.txt"):
            for old_path in output_directory.glob(pattern):
                if old_path.name not in generated:
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
    if not 0 <= args.min_group_margin <= 1:
        raise ValueError("--min-group-margin must be between zero and one")
    if not 0 <= args.min_assigned_fraction <= 1:
        raise ValueError("--min-assigned-fraction must be between zero and one")
    if args.max_rounds < 1:
        raise ValueError("--max-rounds must be positive")
    if args.refinement_rounds < 0:
        raise ValueError("--refinement-rounds must be non-negative")
    if not args.RE:
        raise ValueError("--RE must not be empty")
    if args.min_allelic_block_anchors < 1:
        raise ValueError("--min-allelic-block-anchors must be positive")
    if args.max_allelic_block_configurations < 1:
        raise ValueError("--max-allelic-block-configurations must be positive")

    output_directory = args.output_dir.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    records, order = parse_fasta(args.fasta, args.RE.upper(), args.flank)
    contig_types, inferred_dosage = parse_contig_types(
        args.contig_type, records, args.unknown_dosage_policy
    )
    args.inferred_dosage = inferred_dosage
    dosage = {unitig: DOSAGE_BY_TYPE[contig_types[unitig]] for unitig in records}
    seed_assignments, seed_chromosome = parse_seed_clusters(
        args.clusters_file, records, dosage, args.ploidy
    )
    if args.cluster_assignments:
        cluster_evidence = parse_cluster_evidence(
            args.cluster_assignments, seed_assignments, args.trusted_seed_bases,
            args.min_adjusted_links, args.min_group_margin,
            args.min_assigned_fraction,
        )
        fixed_seed_assignments = {
            unitig: groups for unitig, groups in seed_assignments.items()
            if cluster_evidence[unitig]["trusted"]
        }
    else:
        cluster_evidence = {
            unitig: {"assignment_basis": "not_provided", "trusted": True}
            for unitig in seed_assignments
        }
        fixed_seed_assignments = dict(seed_assignments)
    chromosome = seed_chromosome or args.fasta.name.replace(".putg.fa", "")
    neighbors, raw_link_records, relevant_pairs = load_relevant_links(
        args.full_links, records, dosage, args.hic_link_normalization
    )
    assignments, decisions, round_counts = assign_nonseed_unitigs(
        records, fixed_seed_assignments, seed_assignments, neighbors, dosage, args.ploidy,
        args.min_adjusted_links, args.min_group_margin, args.min_assigned_fraction,
        args.max_rounds, args.low_confidence_policy, args.reviewed_seed_fallback,
    )
    (
        assignments,
        decisions,
        refinement_change_counts,
        refinement_events,
        refinement_stable,
        refinement_oscillation,
    ) = refine_assignments(
        records, assignments, decisions, fixed_seed_assignments, neighbors, dosage,
        args.ploidy, args.min_adjusted_links, args.min_group_margin,
        args.min_assigned_fraction, args.refinement_rounds,
    )
    relaxed_pairs = (
        parse_relaxed_constraints(args.relaxed_constraints)
        if args.relaxed_constraints else set()
    )
    blocks = (
        parse_allelic_blocks(args.allelic_table, records, relaxed_pairs)
        if args.allelic_table else []
    )
    assignments, decisions, block_events = joint_reassign_allelic_blocks(
        records, assignments, decisions, fixed_seed_assignments, neighbors,
        dosage, blocks, args,
    )
    violations = validate(
        records, fixed_seed_assignments, assignments, dosage, block_events
    )
    if violations:
        raise RuntimeError(f"Internal recluster validation failed with {len(violations)} violations")
    summary = write_outputs(
        output_directory, chromosome, records, order, contig_types, dosage,
        seed_assignments, fixed_seed_assignments, cluster_evidence, assignments,
        decisions, round_counts, refinement_change_counts, refinement_events,
        refinement_stable, refinement_oscillation, neighbors, block_events,
        raw_link_records,
        relevant_pairs, violations, args,
    )
    print(
        f"{chromosome}: assigned {summary['assigned']['unitigs']}/{summary['input']['unitigs']} "
        f"unitigs ({summary['assigned']['bp']}/{summary['input']['bp']} bp); "
        f"unassigned={summary['unassigned']['unitigs']}"
    )
    return summary


def parse_arguments():
    parser = argparse.ArgumentParser(description="Dosage-aware PHap v2 chromosome unitig reassignment")
    parser.add_argument("--fasta", required=True, type=Path)
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
    parser.add_argument("--clusters-file", "--clusters_file", required=True, type=Path)
    parser.add_argument("--cluster-assignments", type=Path)
    parser.add_argument(
        "--allelic-table", type=Path,
        help="Chromosome allelic table used for conservative joint block reassignment",
    )
    parser.add_argument(
        "--relaxed-constraints", type=Path,
        help="03.cluster relaxed pairs that must not be reintroduced as strict block edges",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--ploidy", type=int, default=4)
    parser.add_argument("--RE", default="GATC")
    parser.add_argument("--flank", type=int)
    parser.add_argument("--min-adjusted-links", type=float, default=5.0)
    parser.add_argument("--min-group-margin", type=float, default=0.10)
    parser.add_argument("--min-allelic-block-anchors", type=int, default=1)
    parser.add_argument("--max-allelic-block-configurations", type=int, default=256)
    parser.add_argument("--min-assigned-fraction", type=float, default=0.0)
    parser.add_argument("--max-rounds", type=int, default=10)
    parser.add_argument(
        "--refinement-rounds", type=int, default=4,
        help="Synchronous reassignment rounds after initial propagation [4]",
    )
    parser.add_argument("--low-confidence-policy", choices=["defer", "best"], default="defer")
    parser.add_argument("--unknown-dosage-policy", choices=["error", "haplotig"], default="error")
    parser.add_argument(
        "--reviewed-seed-fallback", choices=["retain", "defer"], default="retain",
        help="Retain the original group when Hi-C cannot confidently resolve a reviewed seed [retain]",
    )
    parser.add_argument(
        "--trusted-seed-bases",
        default="hic_supported",
        help="Comma-separated 03.cluster assignment bases retained as immutable anchors",
    )
    args = parser.parse_args()
    args.trusted_seed_bases = {
        basis.strip() for basis in args.trusted_seed_bases.split(",") if basis.strip()
    }
    if not args.trusted_seed_bases:
        parser.error("--trusted-seed-bases must contain at least one assignment basis")
    return args


def main():
    run(parse_arguments())


if __name__ == "__main__":
    main()
