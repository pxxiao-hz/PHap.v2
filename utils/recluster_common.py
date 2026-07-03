"""Shared I/O boundary for chromosome and genome-wide reclustering."""

from __future__ import annotations

import os
import pickle
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Dict, Optional, Tuple

from phap_core.reclustering import ReclusterResult


def read_fasta_with_re_sites(
    fasta_path: str,
    recognition_sequence: str,
) -> Tuple[Dict[str, str], Dict[str, int], Dict[str, int]]:
    if not recognition_sequence:
        raise ValueError("restriction-enzyme recognition sequence must not be empty")
    sequences: Dict[str, list[str]] = {}
    current_id: Optional[str] = None
    with open(fasta_path, encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            if not line.strip():
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                if not current_id or current_id in sequences:
                    raise ValueError(
                        f"{fasta_path}:{line_number}: invalid or duplicate FASTA ID"
                    )
                sequences[current_id] = []
            elif current_id is None:
                raise ValueError(f"{fasta_path}:{line_number}: sequence before FASTA ID")
            else:
                sequences[current_id].append(line.strip())
    joined = {unitig_id: "".join(parts) for unitig_id, parts in sequences.items()}
    return (
        joined,
        {unitig_id: len(sequence) for unitig_id, sequence in joined.items()},
        {
            unitig_id: sequence.upper().count(recognition_sequence.upper())
            for unitig_id, sequence in joined.items()
        },
    )


def load_pair_links(path: str) -> Dict[Tuple[str, str], float]:
    with open(path, "rb") as source:
        value = pickle.load(source)
    if not isinstance(value, dict):
        raise ValueError("Hi-C links pickle must contain a dictionary")
    return value


def parse_group_file(
    path: str,
    *,
    has_count_column: bool,
) -> Tuple[Dict[str, set[str]], Tuple[str, ...]]:
    memberships: Dict[str, set[str]] = {}
    declared_groups = []
    with open(path, encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            minimum = 2 if has_count_column else 1
            if len(fields) < minimum:
                raise ValueError(f"{path}:{line_number}: malformed group row")
            group_id = fields[0]
            if group_id in declared_groups:
                raise ValueError(f"{path}:{line_number}: duplicate group {group_id!r}")
            declared_groups.append(group_id)
            if has_count_column:
                try:
                    expected_count = int(fields[1])
                except ValueError as exc:
                    raise ValueError(
                        f"{path}:{line_number}: invalid group count"
                    ) from exc
                unitigs = fields[2:]
                if len(unitigs) != len(set(unitigs)):
                    raise ValueError(
                        f"{path}:{line_number}: duplicate unitig within group"
                    )
                if expected_count != len(unitigs):
                    raise ValueError(
                        f"{path}:{line_number}: count does not match unitig IDs"
                    )
            else:
                unitigs = fields[1:]
                if len(unitigs) != len(set(unitigs)):
                    raise ValueError(
                        f"{path}:{line_number}: duplicate unitig within group"
                    )
            for unitig_id in unitigs:
                memberships.setdefault(unitig_id, set()).add(group_id)
    if not declared_groups:
        raise ValueError(f"no groups found in {path}")
    return memberships, tuple(declared_groups)


def write_recluster_outputs(
    result: ReclusterResult,
    sequences: Dict[str, str],
    lengths: Dict[str, int],
    restriction_sites: Dict[str, int],
    declared_groups: Iterable[str],
    *,
    output_group_id: Callable[[str], str],
) -> Dict[str, set[str]]:
    memberships = {
        unitig_id: set(group_ids) for unitig_id, group_ids in result.memberships
    }
    group_members: Dict[str, list[str]] = {
        group_id: [] for group_id in sorted(set(declared_groups))
    }
    for unitig_id, group_ids in memberships.items():
        for group_id in group_ids:
            if group_id not in group_members:
                raise ValueError(f"assignment references undeclared group {group_id!r}")
            group_members[group_id].append(unitig_id)

    cluster_lines = []
    combined_fasta: list[str] = []
    mapping_lines = ["output_ID\tsource_ID\tgroup"]
    for group_id in sorted(group_members):
        output_id = output_group_id(group_id)
        unitigs = tuple(sorted(group_members[group_id]))
        cluster_lines.append("\t".join((output_id, *unitigs)))
        _atomic_write_lines(f"{group_id}.reassignment.txt", unitigs)
        fasta_lines: list[str] = []
        metadata_lines = ["#Contig\tRECounts\tLength"]
        for unitig_id in unitigs:
            if unitig_id not in sequences:
                raise ValueError(f"assigned unitig {unitig_id!r} is absent from FASTA")
            fasta_lines.extend((f">{unitig_id}", sequences[unitig_id]))
            metadata_lines.append(
                f"{unitig_id}\t{restriction_sites[unitig_id]}\t{lengths[unitig_id]}"
            )
            renamed = f"{output_id}|{unitig_id}"
            combined_fasta.extend((f">{renamed}", sequences[unitig_id]))
            mapping_lines.append(f"{renamed}\t{unitig_id}\t{output_id}")
        _atomic_write_lines(f"{group_id}.reassignment.fa", fasta_lines)
        _atomic_write_lines(f"{group_id}.txt", metadata_lines)

    _atomic_write_lines("group.reassignment.cluster.txt", cluster_lines)
    _atomic_write_lines("groups.reassignment.fa", combined_fasta)
    _atomic_write_lines("group_sequence_ids.tsv", mapping_lines)
    _write_decision_audits(result)
    unresolved_fasta: list[str] = []
    for decision in result.decisions:
        if decision.status == "assigned":
            continue
        if decision.unitig_id not in sequences:
            raise ValueError(
                f"unresolved unitig {decision.unitig_id!r} is absent from FASTA"
            )
        unresolved_fasta.extend(
            (f">{decision.unitig_id}", sequences[decision.unitig_id])
        )
    _atomic_write_lines("unassigned_unitigs.fa", unresolved_fasta)
    return memberships


def _write_decision_audits(result: ReclusterResult) -> None:
    decision_lines = [
        "unitig_ID\tsource_state\tdosage\tstatus\tselected_groups\t"
        "selected_family\tassigned_locus\treason"
    ]
    score_lines = [
        "unitig_ID\tgroup\traw_score\tre_sites\tre_density\tre_status\t"
        "raw_rank\tre_density_rank"
    ]
    unresolved_lines = ["unitig_ID\tstatus\tassigned_locus\treason"]
    for decision in result.decisions:
        decision_lines.append(
            "\t".join(
                (
                    decision.unitig_id,
                    decision.source_state,
                    "." if decision.dosage is None else str(decision.dosage),
                    decision.status,
                    ",".join(decision.selected_groups),
                    decision.selected_family or ".",
                    decision.assigned_locus or ".",
                    decision.reason,
                )
            )
        )
        if decision.status != "assigned":
            unresolved_lines.append(
                "\t".join(
                    (
                        decision.unitig_id,
                        decision.status,
                        decision.assigned_locus or ".",
                        decision.reason,
                    )
                )
            )
        for score in decision.scores:
            score_lines.append(
                "\t".join(
                    (
                        decision.unitig_id,
                        score.group_id,
                        f"{score.raw_score:.6f}",
                        "." if score.re_sites is None else str(score.re_sites),
                        "." if score.re_density is None else f"{score.re_density:.12g}",
                        score.re_status,
                        str(score.raw_rank),
                        (
                            "."
                            if score.re_density_rank is None
                            else str(score.re_density_rank)
                        ),
                    )
                )
            )
    _atomic_write_lines("recluster_decisions.tsv", decision_lines)
    _atomic_write_lines("recluster_scores.tsv", score_lines)
    _atomic_write_lines("unassigned_unitigs.tsv", unresolved_lines)


def _atomic_write_lines(path: str, lines: Iterable[str]) -> None:
    target = Path(path)
    temporary = target.with_name(f".{target.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as output:
        for line in lines:
            output.write(line)
            output.write("\n")
    os.replace(temporary, target)
