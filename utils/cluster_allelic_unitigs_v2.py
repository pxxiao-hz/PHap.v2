#!/usr/bin/env python3
"""Cluster allelic unitigs under hard dosage and mutual-exclusion constraints."""

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
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from typing import Optional

from dosage import dosage_from_contig_type


class UnsatisfiableConstraintComponent(ValueError):
    def __init__(self, component):
        self.component = tuple(component)
        super().__init__(
            "Allelic constraints are not satisfiable at the configured ploidy"
        )


class ConstraintSearchLimit(RuntimeError):
    def __init__(self, component, max_backtracks):
        self.component = tuple(component)
        super().__init__(
            "Constraint search exceeded "
            f"{max_backtracks} backtracks in one constraint component"
        )


@dataclass(frozen=True)
class TableRow:
    chromosome: str
    start: int
    end: int
    unitigs: tuple[str, ...]


def parse_fasta(path: Path, flank: Optional[int]):
    sequences = {}
    current_id = None
    chunks = []

    def finish_record():
        if current_id is None:
            return
        sequence = "".join(chunks).upper()
        if flank is None:
            re_sequence = sequence
        else:
            flank_size = min(flank, len(sequence) // 2)
            re_sequence = (
                sequence[:flank_size] + sequence[-flank_size:]
                if flank_size > 0
                else ""
            )
        sequences[current_id] = {
            "sequence": sequence,
            "length": len(sequence),
            "re_sites": re_sequence.count("GATC") + 1,
        }

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                finish_record()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                current_id = header.split()[0]
                if current_id in sequences:
                    raise ValueError(f"Duplicate FASTA ID {current_id} in {path}")
                chunks = []
            elif line.strip():
                if current_id is None:
                    raise ValueError(f"Sequence before first FASTA header at {path}:{line_number}")
                chunks.append(line.strip())
    finish_record()
    if not sequences:
        raise ValueError(f"No FASTA records found in {path}")
    return sequences


def parse_contig_types(path: Path):
    types = {}
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
            if len(fields) <= max(id_index, type_index):
                continue
            unitig = fields[id_index]
            if unitig in types:
                raise ValueError(f"Duplicate contig type for {unitig} at line {line_number}")
            types[unitig] = fields[type_index]
    return types


def parse_allelic_table(path: Path, contig_types, ploidy: int):
    rows = []
    chromosome = None
    units = set()
    adjacency = defaultdict(set)
    row_counts = Counter()
    unit_table_bp = Counter()
    edge_overlap_bp = Counter()
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"Malformed allelic table row at {path}:{line_number}")
            row_chromosome = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"Invalid interval at {path}:{line_number}") from exc
            if start < 0 or end <= start:
                raise ValueError(f"Invalid interval at {path}:{line_number}")
            if chromosome is None:
                chromosome = row_chromosome
            elif chromosome != row_chromosome:
                raise ValueError("Per-chromosome cluster input contains multiple chromosomes")
            unitigs = tuple(fields[3:])
            if len(unitigs) != len(set(unitigs)):
                raise ValueError(f"Duplicate unitig within table row at {path}:{line_number}")
            total_dosage = 0
            for unitig in unitigs:
                contig_type = contig_types.get(unitig)
                dosage = dosage_from_contig_type(contig_type)
                if dosage is None:
                    raise ValueError(f"Missing or invalid dosage type for {unitig}: {contig_type}")
                if dosage > ploidy:
                    raise ValueError(f"Dosage for {unitig} exceeds ploidy")
                total_dosage += dosage
                units.add(unitig)
                row_counts[unitig] += 1
                unit_table_bp[unitig] += end - start
            if total_dosage > ploidy:
                raise ValueError(
                    f"Allelic table dosage exceeds ploidy at {path}:{line_number}"
                )
            for unitig1, unitig2 in itertools.combinations(unitigs, 2):
                adjacency[unitig1].add(unitig2)
                adjacency[unitig2].add(unitig1)
                edge_overlap_bp[tuple(sorted((unitig1, unitig2)))] += end - start
            rows.append(TableRow(row_chromosome, start, end, unitigs))
    if not rows:
        raise ValueError(f"No allelic table rows found in {path}")
    for unitig in units:
        adjacency[unitig]
    return (
        chromosome,
        rows,
        units,
        adjacency,
        row_counts,
        unit_table_bp,
        edge_overlap_bp,
    )


def parse_allelic_pair_evidence(path: Optional[Path], chromosome, table_units):
    """Read direct-projection support for table edges.

    The allelic table can contain constraints inferred from a long-path envelope.
    This sidecar preserves which pairs also overlap in accepted alignment blocks.
    """
    evidence = {}
    if path is None:
        return evidence
    required = {
        "target",
        "unitig1",
        "unitig2",
        "direct_projection_overlap_bp",
    }
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"Allelic-pair evidence is missing columns: {', '.join(sorted(missing))}"
            )
        for line_number, row in enumerate(reader, 2):
            if row["target"] != chromosome:
                continue
            unitig1 = row["unitig1"]
            unitig2 = row["unitig2"]
            if unitig1 not in table_units or unitig2 not in table_units:
                continue
            edge = tuple(sorted((unitig1, unitig2)))
            if edge in evidence:
                raise ValueError(
                    f"Duplicate allelic-pair evidence for {edge[0]},{edge[1]} "
                    f"at {path}:{line_number}"
                )
            try:
                direct_overlap_bp = int(row["direct_projection_overlap_bp"])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid direct projection overlap at {path}:{line_number}"
                ) from exc
            if direct_overlap_bp < 0:
                raise ValueError(
                    f"Negative direct projection overlap at {path}:{line_number}"
                )
            evidence[edge] = direct_overlap_bp
    return evidence


def find_overcapacity_cliques(adjacency, dosage, ploidy):
    cliques = []

    def bron_kerbosch(current, candidates, excluded):
        if not candidates and not excluded:
            if sum(dosage[unitig] for unitig in current) > ploidy:
                cliques.append(tuple(sorted(current)))
            return
        pivot_candidates = candidates | excluded
        pivot = max(
            pivot_candidates,
            key=lambda unitig: len(candidates & adjacency[unitig]),
            default=None,
        )
        extension = candidates - (adjacency[pivot] if pivot is not None else set())
        for unitig in sorted(extension):
            neighbors = adjacency[unitig]
            bron_kerbosch(
                current | {unitig},
                candidates & neighbors,
                excluded & neighbors,
            )
            candidates.remove(unitig)
            excluded.add(unitig)

    bron_kerbosch(set(), set(adjacency), set())
    return cliques


def prepare_enforced_constraints(
    original_adjacency,
    dosage,
    ploidy,
    unitig_lengths,
    edge_overlap_bp,
    relaxation_mode,
    protected_edges,
):
    adjacency = {unitig: set(neighbors) for unitig, neighbors in original_adjacency.items()}
    relaxed = []
    while True:
        overcapacity = find_overcapacity_cliques(adjacency, dosage, ploidy)
        if not overcapacity:
            break
        if relaxation_mode == "fail":
            example = overcapacity[0]
            weighted_size = sum(dosage[unitig] for unitig in example)
            raise ValueError(
                "Allelic constraints are not satisfiable at the configured ploidy; "
                f"over-capacity clique dosage={weighted_size}: {', '.join(example)}"
            )
        candidate_edges = set()
        edge_cliques = defaultdict(list)
        for clique in overcapacity:
            for pair in itertools.combinations(clique, 2):
                edge = tuple(sorted(pair))
                candidate_edges.add(edge)
                edge_cliques[edge].append(clique)

        relaxable_edges = candidate_edges - protected_edges
        if not relaxable_edges:
            raise ValueError(
                "Allelic constraints are not satisfiable without removing a protected "
                "high-confidence allelic edge"
            )

        edge = min(
            relaxable_edges,
            key=lambda pair: constraint_edge_weakness(
                pair, unitig_lengths, edge_overlap_bp
            ),
        )
        unitig1, unitig2 = edge
        overlap_bp = edge_overlap_bp[edge]
        relaxed.append(
            {
                "unitig1": unitig1,
                "unitig2": unitig2,
                "overlap_bp": overlap_bp,
                "ratio1": overlap_bp / unitig_lengths[unitig1],
                "ratio2": overlap_bp / unitig_lengths[unitig2],
                "reason": "global_dosage_clique",
                "cliques": [list(clique) for clique in edge_cliques[edge]],
            }
        )
        adjacency[unitig1].remove(unitig2)
        adjacency[unitig2].remove(unitig1)
    return adjacency, relaxed


def constraint_edge_weakness(edge, unitig_lengths, edge_overlap_bp):
    unitig1, unitig2 = edge
    overlap_bp = edge_overlap_bp[edge]
    ratio1 = overlap_bp / unitig_lengths[unitig1]
    ratio2 = overlap_bp / unitig_lengths[unitig2]
    return (min(ratio1, ratio2), overlap_bp, max(ratio1, ratio2), edge)


def relax_weakest_component_edge(
    adjacency,
    component,
    unitig_lengths,
    edge_overlap_bp,
    reason,
    protected_edges,
):
    component = set(component)
    candidate_edges = {
        tuple(sorted((unitig, neighbor)))
        for unitig in component
        for neighbor in adjacency[unitig]
        if neighbor in component
    }
    if not candidate_edges:
        raise ValueError(
            "Constraint component is unsatisfiable but contains no relaxable edges"
        )
    candidate_edges.difference_update(protected_edges)
    if not candidate_edges:
        raise ValueError(
            "Constraint component cannot be solved without removing a protected "
            "high-confidence allelic edge"
        )
    edge = min(
        candidate_edges,
        key=lambda pair: constraint_edge_weakness(
            pair, unitig_lengths, edge_overlap_bp
        ),
    )
    unitig1, unitig2 = edge
    overlap_bp = edge_overlap_bp[edge]
    adjacency[unitig1].remove(unitig2)
    adjacency[unitig2].remove(unitig1)
    return {
        "unitig1": unitig1,
        "unitig2": unitig2,
        "overlap_bp": overlap_bp,
        "ratio1": overlap_bp / unitig_lengths[unitig1],
        "ratio2": overlap_bp / unitig_lengths[unitig2],
        "reason": reason,
        "cliques": [],
    }


def find_protected_long_unitig_edges(
    adjacency,
    lengths,
    edge_overlap_bp,
    min_long_length,
    min_length_ratio,
    min_short_overlap,
):
    protected = set()
    for unitig, neighbors in adjacency.items():
        for neighbor in neighbors:
            edge = tuple(sorted((unitig, neighbor)))
            if edge in protected:
                continue
            length1 = lengths[edge[0]]
            length2 = lengths[edge[1]]
            long_length = max(length1, length2)
            short_length = min(length1, length2)
            if short_length <= 0 or long_length < min_long_length:
                continue
            if long_length / short_length < min_length_ratio:
                continue
            overlap_fraction = edge_overlap_bp[edge] / short_length
            if overlap_fraction >= min_short_overlap:
                protected.add(edge)
    return protected


def find_protected_direct_projection_edges(
    adjacency,
    lengths,
    direct_projection_overlap_bp,
    min_overlap_bp,
    min_short_overlap,
):
    """Protect substantial conflicts supported by accepted alignment blocks."""
    protected = set()
    for edge, overlap_bp in direct_projection_overlap_bp.items():
        unitig1, unitig2 = edge
        if unitig2 not in adjacency.get(unitig1, ()):
            continue
        short_length = min(lengths[unitig1], lengths[unitig2])
        if short_length <= 0 or overlap_bp < min_overlap_bp:
            continue
        if overlap_bp / short_length >= min_short_overlap:
            protected.add(edge)
    return protected


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
    invalid_records = 0
    for key, value in raw_links.items():
        if not isinstance(key, tuple) or len(key) != 2:
            invalid_records += 1
            continue
        unitig1, unitig2 = key
        if unitig1 == unitig2 or unitig1 not in units or unitig2 not in units:
            continue
        try:
            count = float(value)
        except (TypeError, ValueError):
            invalid_records += 1
            continue
        if not math.isfinite(count) or count < 0:
            invalid_records += 1
            continue
        pair = tuple(sorted((unitig1, unitig2)))
        if normalization == "dosage":
            count /= dosage[unitig1] * dosage[unitig2]
        links[pair] += count
    if invalid_records:
        raise ValueError(f"Hi-C link pickle contains {invalid_records} invalid records")
    neighbors = defaultdict(dict)
    for (unitig1, unitig2), count in links.items():
        neighbors[unitig1][unitig2] = count
        neighbors[unitig2][unitig1] = count
    return dict(links), neighbors, len(raw_links)


def masks_for_dosage(ploidy: int, dosage: int):
    masks = []
    for groups in itertools.combinations(range(ploidy), dosage):
        mask = 0
        for group in groups:
            mask |= 1 << group
        masks.append(mask)
    return masks


def groups_in_mask(mask: int, ploidy: int):
    return [group for group in range(ploidy) if mask & (1 << group)]


def mask_size(mask: int):
    return bin(mask).count("1")


class ConstraintClusterer:
    def __init__(
        self,
        units,
        dosage,
        lengths,
        re_sites,
        adjacency,
        links,
        link_neighbors,
        rows,
        ploidy,
        balance_weight,
        max_backtracks,
    ):
        self.units = sorted(units)
        self.dosage = dosage
        self.lengths = lengths
        self.re_sites = re_sites
        self.adjacency = adjacency
        self.links = links
        self.link_neighbors = link_neighbors
        self.rows = rows
        self.ploidy = ploidy
        self.balance_weight = balance_weight
        self.max_backtracks = max_backtracks
        self.domains = {
            unitig: masks_for_dosage(ploidy, dosage[unitig]) for unitig in self.units
        }
        self.backtracks = 0
        self.search_nodes = 0

    def pair_link(self, unitig1, unitig2):
        return self.links.get(tuple(sorted((unitig1, unitig2))), 0.0)

    def select_anchor_row(self):
        candidates = []
        for row_number, row in enumerate(self.rows):
            row_dosage = sum(self.dosage[unitig] for unitig in row.unitigs)
            if row_dosage != self.ploidy:
                continue
            row_set = set(row.unitigs)
            external_links = sum(
                count
                for unitig in row.unitigs
                for neighbor, count in self.link_neighbors.get(unitig, {}).items()
                if neighbor not in row_set
            )
            candidates.append(
                (external_links, row.end - row.start, -row_number, row)
            )
        return max(candidates, default=(None, None, None, None))[3]

    def initial_anchor_assignments(self, anchor_row):
        if anchor_row is None:
            return {}
        assignments = {}
        available = (1 << self.ploidy) - 1
        ordered = sorted(
            anchor_row.unitigs,
            key=lambda unitig: (-self.dosage[unitig], -self.lengths[unitig], unitig),
        )
        for unitig in ordered:
            candidates = [mask for mask in self.domains[unitig] if mask & ~available == 0]
            if not candidates:
                raise ValueError("Unable to initialize the full-dosage anchor row")
            mask = min(candidates)
            assignments[unitig] = mask
            available &= ~mask
        return assignments

    def constraint_components(self):
        components = []
        unseen = set(self.units)
        while unseen:
            seed = min(unseen)
            unseen.remove(seed)
            stack = [seed]
            component = []
            while stack:
                unitig = stack.pop()
                component.append(unitig)
                neighbors = self.adjacency[unitig] & unseen
                unseen.difference_update(neighbors)
                stack.extend(sorted(neighbors, reverse=True))
            components.append(tuple(sorted(component)))
        return components

    def solve(self):
        anchor_row = self.select_anchor_row()
        anchor_assignments = self.initial_anchor_assignments(anchor_row)
        anchor_units = set(anchor_assignments)
        components = self.constraint_components()
        components.sort(
            key=lambda component: (
                not bool(set(component) & anchor_units),
                -sum(len(self.adjacency[unitig]) for unitig in component),
                -sum(self.lengths[unitig] for unitig in component),
                component,
            )
        )

        assignments = {}
        used_anchor = anchor_row
        for component in components:
            fixed = {
                unitig: mask
                for unitig, mask in anchor_assignments.items()
                if unitig in component
            }
            component_assignments = self._solve_component(
                component, fixed, assignments
            )
            if component_assignments is None and fixed:
                component_assignments = self._solve_component(
                    component, {}, assignments
                )
                used_anchor = None
            if component_assignments is None:
                raise UnsatisfiableConstraintComponent(component)
            assignments.update(component_assignments)
        return assignments, used_anchor

    def _solve_component(self, active_units, fixed, base_assignments):
        active_units = tuple(active_units)
        component_start_backtracks = self.backtracks
        assignments = dict(base_assignments)
        assignments.update(fixed)
        group_bp = [0] * self.ploidy
        group_re = [0] * self.ploidy
        for unitig, mask in assignments.items():
            for group in groups_in_mask(mask, self.ploidy):
                group_bp[group] += self.lengths[unitig]
                group_re[group] += self.re_sites[unitig]

        sys.setrecursionlimit(max(2000, len(self.units) * 4))

        def feasible_masks(unitig):
            forbidden = 0
            for neighbor in self.adjacency[unitig]:
                forbidden |= assignments.get(neighbor, 0)
            return [mask for mask in self.domains[unitig] if mask & forbidden == 0]

        def choose_unitig():
            selected = None
            selected_masks = None
            selected_key = None
            for unitig in active_units:
                if unitig in assignments:
                    continue
                masks = feasible_masks(unitig)
                if not masks:
                    return unitig, []
                assigned_neighbor_colors = 0
                for neighbor in self.adjacency[unitig]:
                    assigned_neighbor_colors |= assignments.get(neighbor, 0)
                key = (
                    len(masks),
                    -mask_size(assigned_neighbor_colors),
                    -len(self.adjacency[unitig]),
                    -len(self.link_neighbors.get(unitig, {})),
                    -self.lengths[unitig],
                    unitig,
                )
                if selected_key is None or key < selected_key:
                    selected = unitig
                    selected_masks = masks
                    selected_key = key
            return selected, selected_masks

        def mask_order_key(unitig, mask):
            link_score = 0.0
            for neighbor, count in self.link_neighbors.get(unitig, {}).items():
                neighbor_mask = assignments.get(neighbor, 0)
                link_score += count * mask_size(mask & neighbor_mask)
            projected = list(group_bp)
            for group in groups_in_mask(mask, self.ploidy):
                projected[group] += self.lengths[unitig]
            return (-link_score, max(projected), sum(projected[group] for group in groups_in_mask(mask, self.ploidy)), mask)

        def add(unitig, mask):
            assignments[unitig] = mask
            for group in groups_in_mask(mask, self.ploidy):
                group_bp[group] += self.lengths[unitig]
                group_re[group] += self.re_sites[unitig]

        def remove(unitig, mask):
            del assignments[unitig]
            for group in groups_in_mask(mask, self.ploidy):
                group_bp[group] -= self.lengths[unitig]
                group_re[group] -= self.re_sites[unitig]

        def search():
            if all(unitig in assignments for unitig in active_units):
                return True
            if self.backtracks - component_start_backtracks > self.max_backtracks:
                raise ConstraintSearchLimit(active_units, self.max_backtracks)
            unitig, masks = choose_unitig()
            if not masks:
                self.backtracks += 1
                return False
            self.search_nodes += 1
            for mask in sorted(masks, key=lambda value: mask_order_key(unitig, value)):
                add(unitig, mask)
                forward_ok = all(
                    feasible_masks(neighbor)
                    for neighbor in self.adjacency[unitig]
                    if neighbor not in assignments
                )
                if forward_ok and search():
                    return True
                remove(unitig, mask)
            self.backtracks += 1
            return False

        if not search():
            return None
        return {unitig: assignments[unitig] for unitig in active_units}

    def candidate_score(self, unitig, mask, assignments, group_bp, group_re):
        link_by_group = [0.0] * self.ploidy
        for neighbor, count in self.link_neighbors.get(unitig, {}).items():
            neighbor_mask = assignments.get(neighbor, 0)
            for group in groups_in_mask(neighbor_mask, self.ploidy):
                link_by_group[group] += count
        density = [
            link_by_group[group] / max(group_re[group], 1)
            for group in range(self.ploidy)
        ]
        maximum_density = max(density, default=0.0)
        mean_bp = max(sum(group_bp) / self.ploidy, 1.0)
        score = 0.0
        for group in groups_in_mask(mask, self.ploidy):
            if maximum_density > 0:
                score += density[group] / maximum_density
            projected_load = (group_bp[group] + self.lengths[unitig]) / mean_bp
            score -= self.balance_weight * projected_load
        return score

    def refine(self, assignments, max_rounds):
        group_bp = [0] * self.ploidy
        group_re = [0] * self.ploidy
        for unitig, mask in assignments.items():
            for group in groups_in_mask(mask, self.ploidy):
                group_bp[group] += self.lengths[unitig]
                group_re[group] += self.re_sites[unitig]

        ordered = sorted(
            self.units,
            key=lambda unitig: (
                -sum(self.link_neighbors.get(unitig, {}).values()),
                -len(self.adjacency[unitig]),
                -self.lengths[unitig],
                unitig,
            ),
        )
        moves = 0
        rounds = 0
        for round_number in range(1, max_rounds + 1):
            round_moves = 0
            for unitig in ordered:
                current = assignments[unitig]
                for group in groups_in_mask(current, self.ploidy):
                    group_bp[group] -= self.lengths[unitig]
                    group_re[group] -= self.re_sites[unitig]
                forbidden = 0
                for neighbor in self.adjacency[unitig]:
                    forbidden |= assignments[neighbor]
                candidates = [
                    mask for mask in self.domains[unitig] if mask & forbidden == 0
                ]
                scored = [
                    (self.candidate_score(unitig, mask, assignments, group_bp, group_re), mask)
                    for mask in candidates
                ]
                best_score, best_mask = max(scored, key=lambda item: (item[0], -item[1]))
                current_score = next(score for score, mask in scored if mask == current)
                if best_mask != current and best_score > current_score + 1e-9:
                    assignments[unitig] = best_mask
                    current = best_mask
                    round_moves += 1
                    moves += 1
                for group in groups_in_mask(current, self.ploidy):
                    group_bp[group] += self.lengths[unitig]
                    group_re[group] += self.re_sites[unitig]
            rounds = round_number
            if round_moves == 0:
                break
        return assignments, rounds, moves

    def global_objective(self, assignments):
        group_bp = [0] * self.ploidy
        for unitig, mask in assignments.items():
            for group in groups_in_mask(mask, self.ploidy):
                group_bp[group] += self.lengths[unitig]

        mean_bp = max(sum(group_bp) / self.ploidy, 1.0)
        balance_cv2 = sum(
            ((value - mean_bp) / mean_bp) ** 2 for value in group_bp
        ) / self.ploidy

        supported = 0.0
        possible = 0.0
        for (unitig1, unitig2), count in self.links.items():
            supported += count * mask_size(
                assignments[unitig1] & assignments[unitig2]
            )
            possible += count * min(self.dosage[unitig1], self.dosage[unitig2])
        hic_cohesion = supported / possible if possible else 0.0
        score = hic_cohesion - self.balance_weight * balance_cv2
        return score, hic_cohesion, balance_cv2, group_bp

    def kempe_components(self, assignments, group1, group2):
        bit1 = 1 << group1
        bit2 = 1 << group2
        swappable = {
            unitig
            for unitig, mask in assignments.items()
            if bool(mask & bit1) != bool(mask & bit2)
        }
        seen = set()
        components = []
        for start in sorted(swappable):
            if start in seen:
                continue
            component = []
            stack = [start]
            seen.add(start)
            while stack:
                unitig = stack.pop()
                component.append(unitig)
                for neighbor in self.adjacency[unitig]:
                    if neighbor in swappable and neighbor not in seen:
                        seen.add(neighbor)
                        stack.append(neighbor)
            components.append(tuple(sorted(component)))
        return components

    def refine_kempe(self, assignments, max_rounds):
        initial = self.global_objective(assignments)
        current = initial
        moves = 0
        rounds = 0
        for round_number in range(1, max_rounds + 1):
            best = None
            for group1, group2 in itertools.combinations(range(self.ploidy), 2):
                swap_bits = (1 << group1) | (1 << group2)
                for component in self.kempe_components(assignments, group1, group2):
                    for unitig in component:
                        assignments[unitig] ^= swap_bits
                    candidate = self.global_objective(assignments)
                    for unitig in component:
                        assignments[unitig] ^= swap_bits
                    if candidate[0] <= current[0] + 1e-12:
                        continue
                    key = (candidate[0], candidate[1], -candidate[2])
                    if best is None or key > best[0]:
                        best = (key, group1, group2, component, candidate)
            rounds = round_number
            if best is None:
                break
            _, group1, group2, component, current = best
            swap_bits = (1 << group1) | (1 << group2)
            for unitig in component:
                assignments[unitig] ^= swap_bits
            moves += 1
        return assignments, {
            "rounds": rounds,
            "moves": moves,
            "initial_score": initial[0],
            "final_score": current[0],
            "initial_hic_cohesion": initial[1],
            "final_hic_cohesion": current[1],
            "initial_balance_cv2": initial[2],
            "final_balance_cv2": current[2],
        }

    def edge_overlap_fraction(self, edge, edge_overlap_bp):
        unitig1, unitig2 = edge
        overlap = edge_overlap_bp[edge]
        return max(
            overlap / self.lengths[unitig1],
            overlap / self.lengths[unitig2],
        )

    def relaxable_edges(
        self,
        edges,
        max_relaxations,
        max_overlap,
        edge_overlap_bp,
        protected_edges,
    ):
        if not edges or len(edges) > max_relaxations:
            return False
        return all(
            edge not in protected_edges
            and self.edge_overlap_fraction(edge, edge_overlap_bp) <= max_overlap
            for edge in edges
        )

    def refine_contradicted_constraints(
        self,
        assignments,
        max_rounds,
        min_objective_gain,
        min_adjusted_links,
        min_group_margin,
        max_relaxations,
        max_overlap,
        edge_overlap_bp,
        protected_edges,
    ):
        """Let strong Hi-C evidence override only weak, audited table edges."""
        initial = self.global_objective(assignments)
        current = initial
        moves = []
        relaxed = []
        rounds = 0
        if max_rounds == 0 or max_relaxations == 0 or not self.links:
            return assignments, {
                "rounds": rounds,
                "moves": moves,
                "relaxed_edges": 0,
                "initial_score": initial[0],
                "final_score": initial[0],
                "initial_hic_cohesion": initial[1],
                "final_hic_cohesion": initial[1],
                "initial_balance_cv2": initial[2],
                "final_balance_cv2": initial[2],
            }, relaxed

        assignments = dict(assignments)
        for round_number in range(1, max_rounds + 1):
            best = None
            for unitig in self.units:
                old_mask = assignments[unitig]
                group_re = [0] * self.ploidy
                for other, other_mask in assignments.items():
                    if other == unitig:
                        continue
                    for group in groups_in_mask(other_mask, self.ploidy):
                        group_re[group] += self.re_sites[other]
                group_links = [0.0] * self.ploidy
                for neighbor, count in self.link_neighbors.get(unitig, {}).items():
                    for group in groups_in_mask(
                        assignments[neighbor], self.ploidy
                    ):
                        group_links[group] += count
                densities = [
                    group_links[group] / max(group_re[group], 1)
                    for group in range(self.ploidy)
                ]

                def mask_evidence(mask):
                    selected = groups_in_mask(mask, self.ploidy)
                    unselected = [
                        group
                        for group in range(self.ploidy)
                        if group not in selected
                    ]
                    adjusted_links = sum(group_links[group] for group in selected)
                    if not selected or not unselected:
                        return adjusted_links, 1.0
                    weakest_selected = min(densities[group] for group in selected)
                    strongest_unselected = max(
                        densities[group] for group in unselected
                    )
                    margin = (
                        (weakest_selected - strongest_unselected)
                        / weakest_selected
                        if weakest_selected > 0
                        else 0.0
                    )
                    return adjusted_links, margin

                old_links, old_margin = mask_evidence(old_mask)
                for new_mask in self.domains[unitig]:
                    if new_mask == old_mask:
                        continue
                    violations = tuple(
                        sorted(
                            tuple(sorted((unitig, neighbor)))
                            for neighbor in self.adjacency[unitig]
                            if new_mask & assignments[neighbor]
                        )
                    )
                    if not self.relaxable_edges(
                        violations,
                        max_relaxations,
                        max_overlap,
                        edge_overlap_bp,
                        protected_edges,
                    ):
                        continue
                    new_links, new_margin = mask_evidence(new_mask)
                    if (
                        new_links < min_adjusted_links
                        or new_margin < min_group_margin
                    ):
                        continue

                    assignments[unitig] = new_mask
                    candidate = self.global_objective(assignments)
                    assignments[unitig] = old_mask
                    gain = candidate[0] - current[0]
                    hic_gain = candidate[1] - current[1]
                    if gain < min_objective_gain or hic_gain <= 1e-12:
                        continue
                    key = (
                        candidate[0],
                        candidate[1],
                        -candidate[2],
                        -len(violations),
                        -sum(edge_overlap_bp[edge] for edge in violations),
                        -self.lengths[unitig],
                        unitig,
                        -new_mask,
                    )
                    if best is None or key > best[0]:
                        best = (
                            key,
                            unitig,
                            old_mask,
                            new_mask,
                            violations,
                            candidate,
                            old_links,
                            new_links,
                            old_margin,
                            new_margin,
                        )
            rounds = round_number
            if best is None:
                break

            (
                _,
                unitig,
                old_mask,
                new_mask,
                violations,
                candidate,
                old_links,
                new_links,
                old_margin,
                new_margin,
            ) = best
            assignments[unitig] = new_mask
            move_relaxed = []
            for edge in violations:
                unitig1, unitig2 = edge
                overlap = edge_overlap_bp[edge]
                self.adjacency[unitig1].remove(unitig2)
                self.adjacency[unitig2].remove(unitig1)
                relaxed.append(
                    {
                        "unitig1": unitig1,
                        "unitig2": unitig2,
                        "overlap_bp": overlap,
                        "ratio1": overlap / self.lengths[unitig1],
                        "ratio2": overlap / self.lengths[unitig2],
                        "reason": "hic_contradicted_weak_constraint",
                        "cliques": [],
                    }
                )
                move_relaxed.append(f"{unitig1},{unitig2}")
            moves.append(
                {
                    "round": round_number,
                    "unitig": unitig,
                    "old_groups": ",".join(
                        str(group + 1)
                        for group in groups_in_mask(old_mask, self.ploidy)
                    ),
                    "new_groups": ",".join(
                        str(group + 1)
                        for group in groups_in_mask(new_mask, self.ploidy)
                    ),
                    "relaxed_constraints": ";".join(move_relaxed),
                    "adjusted_links_before": old_links,
                    "adjusted_links_after": new_links,
                    "group_margin_before": old_margin,
                    "group_margin_after": new_margin,
                    "objective_before": current[0],
                    "objective_after": candidate[0],
                    "hic_before": current[1],
                    "hic_after": candidate[1],
                    "balance_cv2_before": current[2],
                    "balance_cv2_after": candidate[2],
                }
            )
            current = candidate

        return assignments, {
            "rounds": rounds,
            "moves": moves,
            "relaxed_edges": len(relaxed),
            "initial_score": initial[0],
            "final_score": current[0],
            "initial_hic_cohesion": initial[1],
            "final_hic_cohesion": current[1],
            "initial_balance_cv2": initial[2],
            "final_balance_cv2": current[2],
        }, relaxed

    def unitig_reference_positions(self):
        weighted_midpoints = defaultdict(float)
        table_bp = Counter()
        for row in self.rows:
            row_bp = row.end - row.start
            midpoint = (row.start + row.end) / 2.0
            for unitig in row.unitigs:
                weighted_midpoints[unitig] += midpoint * row_bp
                table_bp[unitig] += row_bp
        return {
            unitig: weighted_midpoints[unitig] / table_bp[unitig]
            for unitig in self.units
            if table_bp[unitig] > 0
        }

    def refine_phase_blocks(
        self,
        assignments,
        max_rounds,
        min_objective_gain,
        max_boundary_relaxations,
        max_boundary_overlap,
        edge_overlap_bp,
        protected_edges,
    ):
        initial = self.global_objective(assignments)
        current = initial
        moves = []
        relaxed = []
        if max_rounds == 0 or not self.links:
            return assignments, {
                "rounds": 0,
                "moves": moves,
                "relaxed_edges": 0,
                "initial_score": initial[0],
                "final_score": initial[0],
                "initial_hic_cohesion": initial[1],
                "final_hic_cohesion": initial[1],
                "initial_balance_cv2": initial[2],
                "final_balance_cv2": initial[2],
            }, relaxed

        positions = self.unitig_reference_positions()
        position_buckets = defaultdict(list)
        for unitig, position in positions.items():
            position_buckets[position].append(unitig)
        ordered_positions = sorted(position_buckets, reverse=True)
        possible_links = sum(
            count * min(self.dosage[unitig1], self.dosage[unitig2])
            for (unitig1, unitig2), count in self.links.items()
        )
        rounds = 0
        for round_number in range(1, max_rounds + 1):
            best = None
            current_supported = current[1] * possible_links
            for group1, group2 in itertools.combinations(range(self.ploidy), 2):
                bit1 = 1 << group1
                bit2 = 1 << group2
                swap_bits = bit1 | bit2
                suffix = set()
                changed = set()
                violations = set()
                support_delta = 0.0
                group_bp_delta = [0] * self.ploidy

                def swapped_mask(unitig):
                    mask = assignments[unitig]
                    if bool(mask & bit1) != bool(mask & bit2):
                        return mask ^ swap_bits
                    return mask

                for position_number, boundary in enumerate(ordered_positions):
                    for unitig in position_buckets[boundary]:
                        old_mask = assignments[unitig]
                        new_mask = swapped_mask(unitig)
                        for neighbor, count in self.link_neighbors.get(unitig, {}).items():
                            old_shared = mask_size(old_mask & assignments[neighbor])
                            if neighbor in suffix:
                                prior_shared = mask_size(
                                    old_mask & swapped_mask(neighbor)
                                )
                                support_delta -= count * (prior_shared - old_shared)
                            else:
                                new_shared = mask_size(new_mask & assignments[neighbor])
                                support_delta += count * (new_shared - old_shared)
                        for neighbor in self.adjacency[unitig]:
                            edge = tuple(sorted((unitig, neighbor)))
                            if neighbor in suffix:
                                violations.discard(edge)
                            elif new_mask & assignments[neighbor]:
                                violations.add(edge)
                        suffix.add(unitig)
                        if new_mask != old_mask:
                            changed.add(unitig)
                            for group in groups_in_mask(old_mask, self.ploidy):
                                group_bp_delta[group] -= self.lengths[unitig]
                            for group in groups_in_mask(new_mask, self.ploidy):
                                group_bp_delta[group] += self.lengths[unitig]

                    if position_number == len(ordered_positions) - 1 or not changed:
                        continue
                    if len(violations) > max_boundary_relaxations:
                        continue
                    relaxable = True
                    for edge in violations:
                        if edge in protected_edges:
                            relaxable = False
                            break
                        unitig1, unitig2 = edge
                        overlap = edge_overlap_bp[edge]
                        overlap_fraction = max(
                            overlap / self.lengths[unitig1],
                            overlap / self.lengths[unitig2],
                        )
                        if overlap_fraction > max_boundary_overlap:
                            relaxable = False
                            break
                    if not relaxable:
                        continue

                    group_bp = [
                        current[3][group] + group_bp_delta[group]
                        for group in range(self.ploidy)
                    ]
                    mean_bp = max(sum(group_bp) / self.ploidy, 1.0)
                    balance_cv2 = sum(
                        ((value - mean_bp) / mean_bp) ** 2 for value in group_bp
                    ) / self.ploidy
                    hic_cohesion = (
                        (current_supported + support_delta) / possible_links
                        if possible_links
                        else 0.0
                    )
                    score = hic_cohesion - self.balance_weight * balance_cv2
                    candidate_metrics = (
                        score,
                        hic_cohesion,
                        balance_cv2,
                        group_bp,
                    )
                    gain = score - current[0]
                    if gain < min_objective_gain:
                        continue
                    key = (
                        score,
                        hic_cohesion,
                        -balance_cv2,
                        -len(violations),
                        -boundary,
                        -group1,
                        -group2,
                    )
                    if best is None or key > best[0]:
                        best = (
                            key,
                            boundary,
                            group1,
                            group2,
                            tuple(sorted(changed)),
                            tuple(sorted(violations)),
                            candidate_metrics,
                        )
            rounds = round_number
            if best is None:
                break

            (
                _,
                boundary,
                group1,
                group2,
                moved,
                violations,
                candidate_metrics,
            ) = best
            swap_bits = (1 << group1) | (1 << group2)
            assignments = dict(assignments)
            for unitig in moved:
                assignments[unitig] ^= swap_bits
            move_relaxed = []
            for unitig1, unitig2 in violations:
                overlap = edge_overlap_bp[(unitig1, unitig2)]
                self.adjacency[unitig1].remove(unitig2)
                self.adjacency[unitig2].remove(unitig1)
                item = {
                    "unitig1": unitig1,
                    "unitig2": unitig2,
                    "overlap_bp": overlap,
                    "ratio1": overlap / self.lengths[unitig1],
                    "ratio2": overlap / self.lengths[unitig2],
                    "reason": "phase_block_boundary",
                    "cliques": [],
                }
                relaxed.append(item)
                move_relaxed.append(f"{unitig1},{unitig2}")
            moves.append(
                {
                    "round": round_number,
                    "boundary": int(round(boundary)),
                    "group1": group1 + 1,
                    "group2": group2 + 1,
                    "moved_unitigs": len(moved),
                    "relaxed_constraints": ";".join(move_relaxed),
                    "objective_before": current[0],
                    "objective_after": candidate_metrics[0],
                    "hic_before": current[1],
                    "hic_after": candidate_metrics[1],
                    "balance_cv2_before": current[2],
                    "balance_cv2_after": candidate_metrics[2],
                }
            )
            current = candidate_metrics

        return assignments, {
            "rounds": rounds,
            "moves": moves,
            "relaxed_edges": len(relaxed),
            "initial_score": initial[0],
            "final_score": current[0],
            "initial_hic_cohesion": initial[1],
            "final_hic_cohesion": current[1],
            "initial_balance_cv2": initial[2],
            "final_balance_cv2": current[2],
        }, relaxed

    def refine_phase_intervals(
        self,
        assignments,
        max_rounds,
        min_objective_gain,
        max_boundary_relaxations,
        max_boundary_overlap,
        edge_overlap_bp,
        protected_edges,
    ):
        """Swap labels inside two-boundary intervals to escape suffix barriers."""
        initial = self.global_objective(assignments)
        current = initial
        moves = []
        relaxed = []
        if max_rounds == 0 or not self.links:
            return assignments, {
                "rounds": 0,
                "moves": moves,
                "relaxed_edges": 0,
                "initial_score": initial[0],
                "final_score": initial[0],
                "initial_hic_cohesion": initial[1],
                "final_hic_cohesion": initial[1],
                "initial_balance_cv2": initial[2],
                "final_balance_cv2": initial[2],
            }, relaxed

        positions = self.unitig_reference_positions()
        position_buckets = defaultdict(list)
        for unitig, position in positions.items():
            position_buckets[position].append(unitig)
        ordered_positions = sorted(position_buckets)
        possible_links = sum(
            count * min(self.dosage[unitig1], self.dosage[unitig2])
            for (unitig1, unitig2), count in self.links.items()
        )
        assignments = dict(assignments)
        rounds = 0
        for round_number in range(1, max_rounds + 1):
            best = None
            current_supported = current[1] * possible_links
            for group1, group2 in itertools.combinations(range(self.ploidy), 2):
                bit1 = 1 << group1
                bit2 = 1 << group2
                swap_bits = bit1 | bit2

                def swapped_mask(unitig):
                    mask = assignments[unitig]
                    if bool(mask & bit1) != bool(mask & bit2):
                        return mask ^ swap_bits
                    return mask

                active_positions = [
                    position
                    for position in ordered_positions
                    if any(
                        swapped_mask(unitig) != assignments[unitig]
                        for unitig in position_buckets[position]
                    )
                ]

                # Prefixes and suffixes are handled by refine_phase_blocks. Only
                # intervals with sequence on both sides are considered here.
                for left_index in range(1, len(active_positions) - 1):
                    interval = set()
                    changed = set()
                    violations = set()
                    support_delta = 0.0
                    group_bp_delta = [0] * self.ploidy
                    for right_index in range(
                        left_index, len(active_positions) - 1
                    ):
                        right_position = active_positions[right_index]
                        for unitig in position_buckets[right_position]:
                            old_mask = assignments[unitig]
                            new_mask = swapped_mask(unitig)
                            for neighbor, count in self.link_neighbors.get(
                                unitig, {}
                            ).items():
                                old_shared = mask_size(
                                    old_mask & assignments[neighbor]
                                )
                                if neighbor in interval:
                                    prior_shared = mask_size(
                                        old_mask & swapped_mask(neighbor)
                                    )
                                    support_delta -= count * (
                                        prior_shared - old_shared
                                    )
                                else:
                                    new_shared = mask_size(
                                        new_mask & assignments[neighbor]
                                    )
                                    support_delta += count * (
                                        new_shared - old_shared
                                    )
                            for neighbor in self.adjacency[unitig]:
                                edge = tuple(sorted((unitig, neighbor)))
                                if neighbor in interval:
                                    violations.discard(edge)
                                elif new_mask & assignments[neighbor]:
                                    violations.add(edge)
                            interval.add(unitig)
                            if new_mask != old_mask:
                                changed.add(unitig)
                                for group in groups_in_mask(
                                    old_mask, self.ploidy
                                ):
                                    group_bp_delta[group] -= self.lengths[unitig]
                                for group in groups_in_mask(
                                    new_mask, self.ploidy
                                ):
                                    group_bp_delta[group] += self.lengths[unitig]

                        if not changed:
                            continue
                        if violations and not self.relaxable_edges(
                            violations,
                            max_boundary_relaxations,
                            max_boundary_overlap,
                            edge_overlap_bp,
                            protected_edges,
                        ):
                            continue

                        group_bp = [
                            current[3][group] + group_bp_delta[group]
                            for group in range(self.ploidy)
                        ]
                        mean_bp = max(sum(group_bp) / self.ploidy, 1.0)
                        balance_cv2 = sum(
                            ((value - mean_bp) / mean_bp) ** 2
                            for value in group_bp
                        ) / self.ploidy
                        hic_cohesion = (
                            (current_supported + support_delta) / possible_links
                            if possible_links
                            else 0.0
                        )
                        score = hic_cohesion - self.balance_weight * balance_cv2
                        gain = score - current[0]
                        if gain < min_objective_gain:
                            continue
                        candidate = (score, hic_cohesion, balance_cv2, group_bp)
                        left_position = active_positions[left_index]
                        key = (
                            score,
                            hic_cohesion,
                            -balance_cv2,
                            -len(violations),
                            -len(changed),
                            -left_position,
                            -right_position,
                            -group1,
                            -group2,
                        )
                        if best is None or key > best[0]:
                            best = (
                                key,
                                left_position,
                                right_position,
                                group1,
                                group2,
                                tuple(sorted(changed)),
                                tuple(sorted(violations)),
                                candidate,
                            )
            rounds = round_number
            if best is None:
                break

            (
                _,
                interval_start,
                interval_end,
                group1,
                group2,
                moved,
                violations,
                candidate,
            ) = best
            swap_bits = (1 << group1) | (1 << group2)
            for unitig in moved:
                assignments[unitig] ^= swap_bits
            move_relaxed = []
            for unitig1, unitig2 in violations:
                edge = (unitig1, unitig2)
                overlap = edge_overlap_bp[edge]
                self.adjacency[unitig1].remove(unitig2)
                self.adjacency[unitig2].remove(unitig1)
                relaxed.append(
                    {
                        "unitig1": unitig1,
                        "unitig2": unitig2,
                        "overlap_bp": overlap,
                        "ratio1": overlap / self.lengths[unitig1],
                        "ratio2": overlap / self.lengths[unitig2],
                        "reason": "phase_interval_boundary",
                        "cliques": [],
                    }
                )
                move_relaxed.append(f"{unitig1},{unitig2}")
            moves.append(
                {
                    "round": round_number,
                    "interval_start": int(round(interval_start)),
                    "interval_end": int(round(interval_end)),
                    "group1": group1 + 1,
                    "group2": group2 + 1,
                    "moved_unitigs": len(moved),
                    "relaxed_constraints": ";".join(move_relaxed),
                    "objective_before": current[0],
                    "objective_after": candidate[0],
                    "hic_before": current[1],
                    "hic_after": candidate[1],
                    "balance_cv2_before": current[2],
                    "balance_cv2_after": candidate[2],
                }
            )
            current = candidate

        return assignments, {
            "rounds": rounds,
            "moves": moves,
            "relaxed_edges": len(relaxed),
            "initial_score": initial[0],
            "final_score": current[0],
            "initial_hic_cohesion": initial[1],
            "final_hic_cohesion": current[1],
            "initial_balance_cv2": initial[2],
            "final_balance_cv2": current[2],
        }, relaxed


def validate_assignments(assignments, dosage, adjacency, ploidy):
    violations = []
    for unitig, expected_dosage in dosage.items():
        observed = mask_size(assignments.get(unitig, 0))
        if observed != expected_dosage:
            violations.append(
                ("dosage", unitig, "", expected_dosage, observed)
            )
    for unitig1 in sorted(adjacency):
        for unitig2 in sorted(adjacency[unitig1]):
            if unitig1 >= unitig2:
                continue
            overlap = assignments[unitig1] & assignments[unitig2]
            if overlap:
                violations.append(
                    (
                        "allelic_conflict",
                        unitig1,
                        unitig2,
                        "disjoint",
                        ",".join(str(group + 1) for group in groups_in_mask(overlap, ploidy)),
                    )
                )
    return violations


def final_unitig_metrics(
    unitig,
    assignments,
    link_neighbors,
    re_sites,
    ploidy,
):
    group_re = [0] * ploidy
    for other, mask in assignments.items():
        if other == unitig:
            continue
        for group in groups_in_mask(mask, ploidy):
            group_re[group] += re_sites[other]
    group_links = [0.0] * ploidy
    for neighbor, count in link_neighbors.get(unitig, {}).items():
        for group in groups_in_mask(assignments[neighbor], ploidy):
            group_links[group] += count
    densities = [
        group_links[group] / max(group_re[group], 1) for group in range(ploidy)
    ]
    assigned_groups = groups_in_mask(assignments[unitig], ploidy)
    unassigned_groups = [group for group in range(ploidy) if group not in assigned_groups]
    assigned_links = sum(group_links[group] for group in assigned_groups)
    total_links = sum(group_links)
    if unassigned_groups and assigned_groups:
        weakest_group = min(
            assigned_groups, key=lambda group: (densities[group], group)
        )
        strongest_group = max(
            unassigned_groups, key=lambda group: (densities[group], -group)
        )
        weakest_assigned = densities[weakest_group]
        strongest_unassigned = densities[strongest_group]
        margin = (
            (weakest_assigned - strongest_unassigned) / weakest_assigned
            if weakest_assigned > 0
            else 0.0
        )
    else:
        margin = 1.0
    if total_links == 0:
        basis = "constraint_only_no_hic"
    elif assigned_links == 0:
        basis = "constraint_forced_no_assigned_hic"
    elif margin > 0:
        basis = "hic_supported"
    else:
        basis = "constraint_supported_low_margin"
    return group_links, densities, assigned_links, total_links, margin, basis


def write_outputs(
    output_directory,
    chromosome,
    sequences,
    table_units,
    contig_types,
    dosage,
    original_adjacency,
    enforced_adjacency,
    row_counts,
    relaxed_constraints,
    protected_edges,
    protected_edge_reasons,
    direct_projection_overlap_bp,
    assignments,
    links,
    link_neighbors,
    raw_link_records,
    anchor_row,
    clusterer,
    refinement_rounds,
    refinement_moves,
    kempe_metrics,
    constraint_relaxation_metrics,
    post_relax_kempe_metrics,
    phase_block_metrics,
    phase_interval_metrics,
    violations,
    args,
):
    temporary_directory = Path(
        tempfile.mkdtemp(prefix=".cluster.", dir=output_directory)
    )
    try:
        group_units = []
        for group in range(args.ploidy):
            members = sorted(
                unitig
                for unitig in table_units
                if assignments[unitig] & (1 << group)
            )
            group_units.append(members)
            with (temporary_directory / f"g{group + 1}.txt").open("w") as handle:
                for unitig in members:
                    handle.write(unitig + "\n")
            with (temporary_directory / f"g{group + 1}.fa").open("w") as handle:
                for unitig in members:
                    handle.write(f">{unitig}\n{sequences[unitig]['sequence']}\n")

        with (temporary_directory / "group.cluster.txt").open("w") as handle:
            for group, members in enumerate(group_units, 1):
                handle.write(f"group{group}\t{len(members)}\t{' '.join(members)}\n")

        padded_length = max(len(members) for members in group_units)
        combined_stem = "g1g2g3g4" if args.ploidy == 4 else "all_groups"
        with (temporary_directory / f"{combined_stem}.txt").open("w") as handle:
            handle.write("\t".join(f"g{group + 1}" for group in range(args.ploidy)) + "\n")
            for row_number in range(padded_length):
                handle.write(
                    "\t".join(
                        members[row_number] if row_number < len(members) else ""
                        for members in group_units
                    )
                    + "\n"
                )
        with (temporary_directory / f"{combined_stem}.fa").open("w") as handle:
            for group, members in enumerate(group_units, 1):
                for unitig in members:
                    handle.write(f">g{group}_{unitig}\n{sequences[unitig]['sequence']}\n")

        assignment_fields = [
            "unitig",
            "length",
            "contig_type",
            "dosage",
            "groups",
            "group_count",
            "assignment_basis",
            "adjusted_assigned_hic_links",
            "adjusted_total_hic_links",
            "assigned_hic_fraction",
            "hic_density_margin",
            "conflict_degree",
            "table_rows",
            *[f"g{group + 1}_links" for group in range(args.ploidy)],
            *[f"g{group + 1}_density" for group in range(args.ploidy)],
        ]
        basis_counts = Counter()
        with (temporary_directory / "cluster_assignments.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=assignment_fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for unitig in sorted(table_units):
                group_links, densities, assigned_links, total_links, margin, basis = final_unitig_metrics(
                    unitig,
                    assignments,
                    link_neighbors,
                    {item: sequences[item]["re_sites"] for item in table_units},
                    args.ploidy,
                )
                basis_counts[basis] += 1
                groups = [group + 1 for group in groups_in_mask(assignments[unitig], args.ploidy)]
                row = {
                    "unitig": unitig,
                    "length": sequences[unitig]["length"],
                    "contig_type": contig_types[unitig],
                    "dosage": dosage[unitig],
                    "groups": ",".join(map(str, groups)),
                    "group_count": len(groups),
                    "assignment_basis": basis,
                    "adjusted_assigned_hic_links": f"{assigned_links:.6f}",
                    "adjusted_total_hic_links": f"{total_links:.6f}",
                    "assigned_hic_fraction": f"{assigned_links / total_links if total_links else 0.0:.6f}",
                    "hic_density_margin": f"{margin:.6f}",
                    "conflict_degree": len(enforced_adjacency[unitig]),
                    "table_rows": row_counts[unitig],
                }
                row.update(
                    {f"g{group + 1}_links": f"{group_links[group]:.6f}" for group in range(args.ploidy)}
                )
                row.update(
                    {f"g{group + 1}_density": f"{densities[group]:.12f}" for group in range(args.ploidy)}
                )
                writer.writerow(row)

        with (temporary_directory / "cluster_constraint_violations.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["violation", "unitig1", "unitig2", "expected", "observed"])
            writer.writerows(violations)

        relaxed_fields = [
            "unitig1",
            "unitig2",
            "overlap_bp",
            "ratio1",
            "ratio2",
            "reason",
            "source_cliques",
            "final_shared_groups",
        ]
        relaxed_shared = 0
        with (temporary_directory / "cluster_relaxed_constraints.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=relaxed_fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for item in relaxed_constraints:
                shared = groups_in_mask(
                    assignments[item["unitig1"]] & assignments[item["unitig2"]],
                    args.ploidy,
                )
                if shared:
                    relaxed_shared += 1
                writer.writerow(
                    {
                        "unitig1": item["unitig1"],
                        "unitig2": item["unitig2"],
                        "overlap_bp": item["overlap_bp"],
                        "ratio1": f'{item["ratio1"]:.6f}',
                        "ratio2": f'{item["ratio2"]:.6f}',
                        "reason": item["reason"],
                        "source_cliques": ";".join(
                            ",".join(clique) for clique in item["cliques"]
                        ),
                        "final_shared_groups": ",".join(str(group + 1) for group in shared),
                    }
                )

        with (temporary_directory / "cluster_protected_constraints.tsv").open(
            "w", newline=""
        ) as handle:
            fields = [
                "unitig1",
                "unitig2",
                "protection_reason",
                "long_unitig",
                "short_unitig",
                "long_length",
                "short_length",
                "overlap_bp",
                "short_overlap_fraction",
                "direct_projection_overlap_bp",
                "direct_projection_short_fraction",
            ]
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            protected_overlap = Counter()
            for row in clusterer.rows:
                for pair in itertools.combinations(row.unitigs, 2):
                    edge = tuple(sorted(pair))
                    if edge in protected_edges:
                        protected_overlap[edge] += row.end - row.start
            for unitig1, unitig2 in sorted(protected_edges):
                long_unitig, short_unitig = sorted(
                    (unitig1, unitig2),
                    key=lambda unitig: (sequences[unitig]["length"], unitig),
                    reverse=True,
                )
                long_length = sequences[long_unitig]["length"]
                short_length = sequences[short_unitig]["length"]
                overlap_bp = protected_overlap[(unitig1, unitig2)]
                direct_overlap_bp = direct_projection_overlap_bp.get(
                    (unitig1, unitig2), 0
                )
                writer.writerow(
                    {
                        "unitig1": unitig1,
                        "unitig2": unitig2,
                        "protection_reason": ",".join(
                            sorted(protected_edge_reasons[(unitig1, unitig2)])
                        ),
                        "long_unitig": long_unitig,
                        "short_unitig": short_unitig,
                        "long_length": long_length,
                        "short_length": short_length,
                        "overlap_bp": overlap_bp,
                        "short_overlap_fraction": f"{overlap_bp / short_length:.6f}",
                        "direct_projection_overlap_bp": direct_overlap_bp,
                        "direct_projection_short_fraction": (
                            f"{direct_overlap_bp / short_length:.6f}"
                        ),
                    }
                )

        phase_fields = [
            "round",
            "boundary",
            "group1",
            "group2",
            "moved_unitigs",
            "relaxed_constraints",
            "objective_before",
            "objective_after",
            "hic_before",
            "hic_after",
            "balance_cv2_before",
            "balance_cv2_after",
        ]
        with (temporary_directory / "cluster_phase_block_moves.tsv").open(
            "w", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=phase_fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for move in phase_block_metrics["moves"]:
                row = dict(move)
                for field in (
                    "objective_before",
                    "objective_after",
                    "hic_before",
                    "hic_after",
                    "balance_cv2_before",
                    "balance_cv2_after",
                ):
                    row[field] = f'{row[field]:.12f}'
                writer.writerow(row)

        constraint_move_fields = [
            "round",
            "unitig",
            "old_groups",
            "new_groups",
            "relaxed_constraints",
            "adjusted_links_before",
            "adjusted_links_after",
            "group_margin_before",
            "group_margin_after",
            "objective_before",
            "objective_after",
            "hic_before",
            "hic_after",
            "balance_cv2_before",
            "balance_cv2_after",
        ]
        with (temporary_directory / "cluster_constraint_relaxation_moves.tsv").open(
            "w", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=constraint_move_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for move in constraint_relaxation_metrics["moves"]:
                row = dict(move)
                for field in (
                    "adjusted_links_before",
                    "adjusted_links_after",
                    "group_margin_before",
                    "group_margin_after",
                    "objective_before",
                    "objective_after",
                    "hic_before",
                    "hic_after",
                    "balance_cv2_before",
                    "balance_cv2_after",
                ):
                    row[field] = f'{row[field]:.12f}'
                writer.writerow(row)

        interval_fields = [
            "round",
            "interval_start",
            "interval_end",
            "group1",
            "group2",
            "moved_unitigs",
            "relaxed_constraints",
            "objective_before",
            "objective_after",
            "hic_before",
            "hic_after",
            "balance_cv2_before",
            "balance_cv2_after",
        ]
        with (temporary_directory / "cluster_phase_interval_moves.tsv").open(
            "w", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=interval_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for move in phase_interval_metrics["moves"]:
                row = dict(move)
                for field in (
                    "objective_before",
                    "objective_after",
                    "hic_before",
                    "hic_after",
                    "balance_cv2_before",
                    "balance_cv2_after",
                ):
                    row[field] = f'{row[field]:.12f}'
                writer.writerow(row)

        group_summary = {}
        for group, members in enumerate(group_units, 1):
            group_summary[f"g{group}"] = {
                "unitigs": len(members),
                "bp": sum(sequences[unitig]["length"] for unitig in members),
                "re_sites": sum(sequences[unitig]["re_sites"] for unitig in members),
            }
        summary = {
            "chromosome": chromosome,
            "inputs": {
                "fasta": str(args.fasta.resolve()),
                "full_links": str(args.full_links.resolve()),
                "contig_type": (
                    str(args.contig_type.resolve()) if args.contig_type else None
                ),
                "allelic_table": str(args.allelic_table.resolve()),
                "allelic_pairs": (
                    str(args.allelic_pairs.resolve()) if args.allelic_pairs else None
                ),
            },
            "parameters": {
                "ploidy": args.ploidy,
                "no_collapse": args.no_collapse,
                "include_all_fasta_unitigs": args.include_all_fasta_unitigs,
                "flank": args.flank,
                "hic_link_normalization": args.hic_link_normalization,
                "balance_weight": args.balance_weight,
                "max_refinement_rounds": args.max_refinement_rounds,
                "max_phase_block_rounds": args.max_phase_block_rounds,
                "min_phase_block_gain": args.min_phase_block_gain,
                "max_phase_boundary_relaxations": args.max_phase_boundary_relaxations,
                "max_phase_boundary_overlap": args.max_phase_boundary_overlap,
                "max_constraint_relaxation_rounds": args.max_constraint_relaxation_rounds,
                "min_constraint_relaxation_gain": args.min_constraint_relaxation_gain,
                "min_constraint_relaxation_links": args.min_constraint_relaxation_links,
                "min_constraint_relaxation_margin": args.min_constraint_relaxation_margin,
                "max_constraint_relaxations_per_move": args.max_constraint_relaxations_per_move,
                "max_constraint_relaxation_overlap": args.max_constraint_relaxation_overlap,
                "max_phase_interval_rounds": args.max_phase_interval_rounds,
                "max_phase_interval_relaxations": args.max_phase_interval_relaxations,
                "max_backtracks": args.max_backtracks,
                "constraint_relaxation": args.constraint_relaxation,
                "protected_long_unitig_length": args.protected_long_unitig_length,
                "protected_length_ratio": args.protected_length_ratio,
                "protected_short_overlap": args.protected_short_overlap,
                "protected_direct_overlap_bp": args.protected_direct_overlap_bp,
                "protected_direct_short_overlap": args.protected_direct_short_overlap,
            },
            "table": {
                "rows": len(clusterer.rows),
                "unitigs": len(table_units),
                "row_unitigs": len(
                    {unitig for row in clusterer.rows for unitig in row.unitigs}
                ),
                "original_conflict_pairs": sum(len(value) for value in original_adjacency.values()) // 2,
                "enforced_conflict_pairs": sum(len(value) for value in enforced_adjacency.values()) // 2,
                "relaxed_conflict_pairs": len(relaxed_constraints),
                "protected_conflict_pairs": len(protected_edges),
                "protected_direct_projection_pairs": sum(
                    "direct_projection" in reasons
                    for reasons in protected_edge_reasons.values()
                ),
            },
            "hic": {
                "raw_link_records": raw_link_records,
                "relevant_pairs": len(links),
            },
            "constraint_search": {
                "search_nodes": clusterer.search_nodes,
                "backtracks": clusterer.backtracks,
                "anchor": (
                    {
                        "start": anchor_row.start,
                        "end": anchor_row.end,
                        "unitigs": list(anchor_row.unitigs),
                    }
                    if anchor_row is not None
                    else None
                ),
            },
            "refinement": {
                "rounds": (
                    refinement_rounds
                    + kempe_metrics["rounds"]
                    + constraint_relaxation_metrics["rounds"]
                    + post_relax_kempe_metrics["rounds"]
                    + phase_block_metrics["rounds"]
                    + phase_interval_metrics["rounds"]
                ),
                "moves": (
                    refinement_moves
                    + kempe_metrics["moves"]
                    + len(constraint_relaxation_metrics["moves"])
                    + post_relax_kempe_metrics["moves"]
                    + len(phase_block_metrics["moves"])
                    + len(phase_interval_metrics["moves"])
                ),
                "single_unit_rounds": refinement_rounds,
                "single_unit_moves": refinement_moves,
                "kempe_rounds": kempe_metrics["rounds"],
                "kempe_moves": kempe_metrics["moves"],
                "constraint_relaxation_rounds": constraint_relaxation_metrics["rounds"],
                "constraint_relaxation_moves": len(constraint_relaxation_metrics["moves"]),
                "constraint_relaxation_relaxed_edges": constraint_relaxation_metrics["relaxed_edges"],
                "post_relax_kempe_rounds": post_relax_kempe_metrics["rounds"],
                "post_relax_kempe_moves": post_relax_kempe_metrics["moves"],
                "objective_initial": kempe_metrics["initial_score"],
                "objective_after_kempe": kempe_metrics["final_score"],
                "objective_after_constraint_relaxation": constraint_relaxation_metrics["final_score"],
                "objective_after_post_relax_kempe": post_relax_kempe_metrics["final_score"],
                "objective_after_phase_blocks": phase_block_metrics["final_score"],
                "objective_final": phase_interval_metrics["final_score"],
                "hic_cohesion_initial": kempe_metrics["initial_hic_cohesion"],
                "hic_cohesion_after_kempe": kempe_metrics["final_hic_cohesion"],
                "hic_cohesion_final": phase_interval_metrics["final_hic_cohesion"],
                "balance_cv2_initial": kempe_metrics["initial_balance_cv2"],
                "balance_cv2_after_kempe": kempe_metrics["final_balance_cv2"],
                "balance_cv2_final": phase_interval_metrics["final_balance_cv2"],
                "phase_block_rounds": phase_block_metrics["rounds"],
                "phase_block_moves": len(phase_block_metrics["moves"]),
                "phase_block_relaxed_edges": phase_block_metrics["relaxed_edges"],
                "phase_block_objective_initial": phase_block_metrics["initial_score"],
                "phase_block_objective_final": phase_block_metrics["final_score"],
                "phase_interval_rounds": phase_interval_metrics["rounds"],
                "phase_interval_moves": len(phase_interval_metrics["moves"]),
                "phase_interval_relaxed_edges": phase_interval_metrics["relaxed_edges"],
                "phase_interval_objective_initial": phase_interval_metrics["initial_score"],
                "phase_interval_objective_final": phase_interval_metrics["final_score"],
            },
            "assignment_basis_counts": dict(sorted(basis_counts.items())),
            "groups": group_summary,
            "validation": {
                "allelic_conflicts": sum(item[0] == "allelic_conflict" for item in violations),
                "dosage_errors": sum(item[0] == "dosage" for item in violations),
                "unassigned_table_unitigs": sum(unitig not in assignments for unitig in table_units),
                "relaxed_constraints_sharing_final_group": relaxed_shared,
            },
        }
        with (temporary_directory / "cluster_summary.json").open("w") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")

        generated = {path.name for path in temporary_directory.iterdir()}
        for old_name in (
            "full.links.txt",
            "flank.links.txt",
            "all_groups.txt" if args.ploidy == 4 else "g1g2g3g4.txt",
            "all_groups.fa" if args.ploidy == 4 else "g1g2g3g4.fa",
        ):
            old_path = output_directory / old_name
            if old_path.exists():
                old_path.unlink()
        for old_path in output_directory.glob("g*.*"):
            if re.fullmatch(r"g[1-9][0-9]*\.(?:fa|txt)", old_path.name) and old_path.name not in generated:
                old_path.unlink()
        for path in temporary_directory.iterdir():
            os.replace(path, output_directory / path.name)
        return summary
    finally:
        shutil.rmtree(temporary_directory, ignore_errors=True)


def run(args):
    if args.ploidy < 2:
        raise ValueError("--ploidy must be at least two")
    if args.flank is not None and args.flank < 0:
        raise ValueError("--flank must be non-negative")
    if args.include_all_fasta_unitigs and not args.no_collapse:
        raise ValueError("--include-all-fasta-unitigs currently requires --no-collapse")
    if args.balance_weight < 0:
        raise ValueError("--balance-weight must be non-negative")
    if args.max_refinement_rounds < 0:
        raise ValueError("--max-refinement-rounds must be non-negative")
    if args.max_phase_block_rounds < 0:
        raise ValueError("--max-phase-block-rounds must be non-negative")
    if args.min_phase_block_gain < 0:
        raise ValueError("--min-phase-block-gain must be non-negative")
    if args.max_phase_boundary_relaxations < 0:
        raise ValueError("--max-phase-boundary-relaxations must be non-negative")
    if not 0 <= args.max_phase_boundary_overlap <= 1:
        raise ValueError("--max-phase-boundary-overlap must be between zero and one")
    if args.max_constraint_relaxation_rounds < 0:
        raise ValueError("--max-constraint-relaxation-rounds must be non-negative")
    if args.min_constraint_relaxation_gain < 0:
        raise ValueError("--min-constraint-relaxation-gain must be non-negative")
    if args.min_constraint_relaxation_links < 0:
        raise ValueError("--min-constraint-relaxation-links must be non-negative")
    if not 0 <= args.min_constraint_relaxation_margin <= 1:
        raise ValueError(
            "--min-constraint-relaxation-margin must be between zero and one"
        )
    if args.max_constraint_relaxations_per_move < 0:
        raise ValueError("--max-constraint-relaxations-per-move must be non-negative")
    if not 0 <= args.max_constraint_relaxation_overlap <= 1:
        raise ValueError(
            "--max-constraint-relaxation-overlap must be between zero and one"
        )
    if args.max_phase_interval_rounds < 0:
        raise ValueError("--max-phase-interval-rounds must be non-negative")
    if args.max_phase_interval_relaxations < 0:
        raise ValueError("--max-phase-interval-relaxations must be non-negative")
    if args.protected_long_unitig_length < 0:
        raise ValueError("--protected-long-unitig-length must be non-negative")
    if args.protected_length_ratio < 1:
        raise ValueError("--protected-length-ratio must be at least one")
    if not 0 <= args.protected_short_overlap <= 1:
        raise ValueError("--protected-short-overlap must be between zero and one")
    if args.protected_direct_overlap_bp < 0:
        raise ValueError("--protected-direct-overlap-bp must be non-negative")
    if not 0 <= args.protected_direct_short_overlap <= 1:
        raise ValueError(
            "--protected-direct-short-overlap must be between zero and one"
        )
    output_directory = args.output_dir.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    sequences = parse_fasta(args.fasta, args.flank)
    if args.no_collapse and args.contig_type is not None:
        raise ValueError("--no-collapse and --contig-type are mutually exclusive")
    if not args.no_collapse and args.contig_type is None:
        raise ValueError("--contig-type is required unless --no-collapse is used")
    contig_types = (
        {unitig: "haplotig" for unitig in sequences}
        if args.no_collapse
        else parse_contig_types(args.contig_type)
    )
    (
        chromosome,
        rows,
        table_units,
        original_adjacency,
        row_counts,
        unit_table_bp,
        edge_overlap_bp,
    ) = parse_allelic_table(args.allelic_table, contig_types, args.ploidy)
    missing_fasta = sorted(table_units - set(sequences))
    if missing_fasta:
        raise ValueError(
            f"{len(missing_fasta)} table unitigs are missing from chromosome FASTA: "
            + ", ".join(missing_fasta[:5])
        )
    if args.include_all_fasta_unitigs:
        table_units.update(sequences)
        for unitig in table_units:
            original_adjacency[unitig]
            row_counts[unitig]
            unit_table_bp[unitig]
    dosage = {
        unitig: dosage_from_contig_type(contig_types[unitig])
        for unitig in table_units
    }
    lengths = {unitig: sequences[unitig]["length"] for unitig in table_units}
    direct_projection_overlap_bp = parse_allelic_pair_evidence(
        args.allelic_pairs, chromosome, table_units
    )
    protected_long_edges = find_protected_long_unitig_edges(
        original_adjacency,
        lengths,
        edge_overlap_bp,
        args.protected_long_unitig_length,
        args.protected_length_ratio,
        args.protected_short_overlap,
    )
    protected_direct_edges = find_protected_direct_projection_edges(
        original_adjacency,
        lengths,
        direct_projection_overlap_bp,
        args.protected_direct_overlap_bp,
        args.protected_direct_short_overlap,
    )
    protected_edges = protected_long_edges | protected_direct_edges
    protected_edge_reasons = defaultdict(set)
    for edge in protected_long_edges:
        protected_edge_reasons[edge].add("long_unitig_containment")
    for edge in protected_direct_edges:
        protected_edge_reasons[edge].add("direct_projection")
    adjacency, relaxed_constraints = prepare_enforced_constraints(
        original_adjacency,
        dosage,
        args.ploidy,
        lengths,
        edge_overlap_bp,
        args.constraint_relaxation,
        protected_edges,
    )
    re_sites = {unitig: sequences[unitig]["re_sites"] for unitig in table_units}
    links, link_neighbors, raw_link_records = load_relevant_links(
        args.full_links, table_units, dosage, args.hic_link_normalization
    )
    failed_search_nodes = 0
    failed_backtracks = 0
    while True:
        clusterer = ConstraintClusterer(
            table_units,
            dosage,
            lengths,
            re_sites,
            adjacency,
            links,
            link_neighbors,
            rows,
            args.ploidy,
            args.balance_weight,
            args.max_backtracks,
        )
        try:
            assignments, anchor_row = clusterer.solve()
            break
        except (UnsatisfiableConstraintComponent, ConstraintSearchLimit) as exc:
            failed_search_nodes += clusterer.search_nodes
            failed_backtracks += clusterer.backtracks
            if args.constraint_relaxation == "fail":
                raise ValueError(
                    "Allelic constraints are not satisfiable at the configured ploidy"
                ) from exc
            reason = (
                "constraint_search_limit"
                if isinstance(exc, ConstraintSearchLimit)
                else "constraint_component_unsatisfiable"
            )
            relaxed_constraints.append(
                relax_weakest_component_edge(
                    adjacency,
                    exc.component,
                    lengths,
                    edge_overlap_bp,
                    reason,
                    protected_edges,
                )
            )
    clusterer.search_nodes += failed_search_nodes
    clusterer.backtracks += failed_backtracks
    assignments, refinement_rounds, refinement_moves = clusterer.refine(
        assignments, args.max_refinement_rounds
    )
    assignments, kempe_metrics = clusterer.refine_kempe(
        assignments, args.max_refinement_rounds
    )
    assignments, constraint_relaxation_metrics, contradicted_relaxed = (
        clusterer.refine_contradicted_constraints(
            assignments,
            args.max_constraint_relaxation_rounds,
            args.min_constraint_relaxation_gain,
            args.min_constraint_relaxation_links,
            args.min_constraint_relaxation_margin,
            args.max_constraint_relaxations_per_move,
            args.max_constraint_relaxation_overlap,
            edge_overlap_bp,
            protected_edges,
        )
    )
    relaxed_constraints.extend(contradicted_relaxed)
    assignments, post_relax_kempe_metrics = clusterer.refine_kempe(
        assignments, args.max_refinement_rounds
    )
    assignments, phase_block_metrics, phase_relaxed = clusterer.refine_phase_blocks(
        assignments,
        args.max_phase_block_rounds,
        args.min_phase_block_gain,
        args.max_phase_boundary_relaxations,
        args.max_phase_boundary_overlap,
        edge_overlap_bp,
        protected_edges,
    )
    relaxed_constraints.extend(phase_relaxed)
    assignments, phase_interval_metrics, interval_relaxed = (
        clusterer.refine_phase_intervals(
            assignments,
            args.max_phase_interval_rounds,
            args.min_phase_block_gain,
            args.max_phase_interval_relaxations,
            args.max_phase_boundary_overlap,
            edge_overlap_bp,
            protected_edges,
        )
    )
    relaxed_constraints.extend(interval_relaxed)
    violations = validate_assignments(assignments, dosage, adjacency, args.ploidy)
    if violations:
        raise RuntimeError(
            f"Internal cluster validation failed with {len(violations)} hard-constraint violations"
        )
    summary = write_outputs(
        output_directory,
        chromosome,
        sequences,
        table_units,
        contig_types,
        dosage,
        original_adjacency,
        adjacency,
        row_counts,
        relaxed_constraints,
        protected_edges,
        protected_edge_reasons,
        direct_projection_overlap_bp,
        assignments,
        links,
        link_neighbors,
        raw_link_records,
        anchor_row,
        clusterer,
        refinement_rounds,
        refinement_moves,
        kempe_metrics,
        constraint_relaxation_metrics,
        post_relax_kempe_metrics,
        phase_block_metrics,
        phase_interval_metrics,
        violations,
        args,
    )
    counts = ", ".join(
        f'{group}={metrics["unitigs"]}' for group, metrics in summary["groups"].items()
    )
    print(f"{chromosome}: clustered {len(table_units)} table unitigs; {counts}")
    return summary


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Constraint-safe PHap v2 allelic-unitig clustering"
    )
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--full-links", required=True, type=Path)
    parser.add_argument(
        "--hic-link-normalization",
        choices=["dosage", "raw"],
        default="dosage",
        help=(
            "Normalize each Hi-C count by the product of the two unitig dosages, "
            "or use raw read-pair counts [dosage]"
        ),
    )
    parser.add_argument("--contig-type", type=Path)
    parser.add_argument(
        "--no-collapse",
        action="store_true",
        help="Treat every unitig as single-copy dosage one; --contig-type is not required",
    )
    parser.add_argument(
        "--include-all-fasta-unitigs",
        action="store_true",
        help=(
            "Cluster every chromosome-FASTA unitig, including unitigs absent "
            "from allelic-table rows"
        ),
    )
    parser.add_argument("--allelic-table", required=True, type=Path)
    parser.add_argument(
        "--allelic-pairs",
        type=Path,
        help=(
            "allelic_pairs.tsv sidecar containing direct-projection evidence; "
            "enables protection of high-confidence same-locus edges"
        ),
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--ploidy", type=int, default=4)
    parser.add_argument("--flank", type=int)
    parser.add_argument("--balance-weight", type=float, default=1.0)
    parser.add_argument("--max-refinement-rounds", type=int, default=10)
    parser.add_argument("--max-phase-block-rounds", type=int, default=10)
    parser.add_argument("--min-phase-block-gain", type=float, default=0.005)
    parser.add_argument("--max-phase-boundary-relaxations", type=int, default=1)
    parser.add_argument("--max-phase-boundary-overlap", type=float, default=0.30)
    parser.add_argument("--max-constraint-relaxation-rounds", type=int, default=50)
    parser.add_argument("--min-constraint-relaxation-gain", type=float, default=0.0)
    parser.add_argument("--min-constraint-relaxation-links", type=float, default=5.0)
    parser.add_argument("--min-constraint-relaxation-margin", type=float, default=0.10)
    parser.add_argument("--max-constraint-relaxations-per-move", type=int, default=2)
    parser.add_argument("--max-constraint-relaxation-overlap", type=float, default=0.30)
    parser.add_argument("--max-phase-interval-rounds", type=int, default=4)
    parser.add_argument("--max-phase-interval-relaxations", type=int, default=2)
    parser.add_argument("--max-backtracks", type=int, default=1_000_000)
    parser.add_argument("--protected-long-unitig-length", type=int, default=5000000)
    parser.add_argument("--protected-length-ratio", type=float, default=5.0)
    parser.add_argument("--protected-short-overlap", type=float, default=0.50)
    parser.add_argument("--protected-direct-overlap-bp", type=int, default=1_000_000)
    parser.add_argument(
        "--protected-direct-short-overlap", type=float, default=0.20
    )
    parser.add_argument(
        "--constraint-relaxation",
        choices=["weakest", "fail"],
        default="weakest",
        help="Resolve global over-ploidy cliques by removing their weakest edge, or fail",
    )
    return parser.parse_args()


def main():
    run(parse_arguments())


if __name__ == "__main__":
    main()
