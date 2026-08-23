#!/usr/bin/env python3
"""Build an allelic unitig table from actual reference-coordinate overlap.

The v1 implementation grouped every unitig touching the same fixed-size bin.
That makes adjacent, non-overlapping unitigs mutually exclusive. This version
uses a sweep over merged alignment intervals, so unitigs share a row only when
their projected reference intervals overlap by a real, configurable amount.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from itertools import combinations
from pathlib import Path


KNOWN_DOSAGE = {
    "haplotig": 1,
    "diplotig": 2,
    "triplotig": 3,
    "tetraplotig": 4,
}


@dataclass(frozen=True)
class Projection:
    target: str
    unitig: str
    start: int
    end: int
    aligned_bp: int
    projection_coverage: float
    identity: float
    query_coverage: float
    query_span_coverage: float
    confidence: float
    dosage: int
    contig_type: str
    placement_mode: str
    block_index: int
    block_count: int
    target_span_bp: int
    query_aligned_bp: int
    query_span_bp: int
    target_aligned_fraction: float
    query_aligned_fraction: float
    chain_identity: float
    mean_mapq: float
    strand: str
    collinearity: float
    projection_class: str
    query_length: int
    constraint_role: str


@dataclass(frozen=True)
class Alignment:
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_start: int
    target_end: int
    matches: int
    block_length: int
    mapq: int


@dataclass(frozen=True)
class Segment:
    target: str
    start: int
    end: int
    unitigs: tuple[str, ...]


def union_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


def interval_length(intervals) -> int:
    return sum(end - start for start, end in intervals)


def intersection_length(intervals1, intervals2) -> int:
    """Return union-aware overlap between two interval collections."""
    left = union_intervals(intervals1)
    right = union_intervals(intervals2)
    left_index = 0
    right_index = 0
    overlap_bp = 0
    while left_index < len(left) and right_index < len(right):
        start = max(left[left_index][0], right[right_index][0])
        end = min(left[left_index][1], right[right_index][1])
        overlap_bp += max(0, end - start)
        if left[left_index][1] < right[right_index][1]:
            left_index += 1
        else:
            right_index += 1
    return overlap_bp


def parse_contig_types(path: Path):
    contig_types = {}
    with path.open() as handle:
        header = handle.readline().split()
        try:
            id_index = header.index("contig_ID")
            type_index = header.index("contig_type")
        except ValueError as exc:
            raise ValueError(
                "contig type file must contain contig_ID and contig_type columns"
            ) from exc
        for line in handle:
            fields = line.split()
            if len(fields) > max(id_index, type_index):
                contig_types[fields[id_index]] = fields[type_index]
    return contig_types


def parse_gfa_links(path: Path | None):
    links = set()
    if path is None:
        return links
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.startswith("L\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"Malformed GFA L record at line {line_number}")
            links.add(tuple(sorted((fields[1], fields[3]))))
    return links


def parse_chain_paf(path: Path):
    records = defaultdict(list)
    query_lengths = {}
    matches = Counter()
    block_lengths = Counter()
    query_intervals = defaultdict(list)
    target_intervals = defaultdict(list)
    placement_modes = defaultdict(set)

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF line {line_number}")
            query = fields[0]
            target = fields[5]
            key = (target, query)
            query_lengths[query] = int(fields[1])
            matches[key] += int(fields[9])
            block_lengths[key] += int(fields[10])
            query_intervals[key].append((int(fields[2]), int(fields[3])))
            target_interval = (int(fields[7]), int(fields[8]))
            target_intervals[key].append(target_interval)
            records[key].append(
                Alignment(
                    query_length=int(fields[1]),
                    query_start=int(fields[2]),
                    query_end=int(fields[3]),
                    strand=fields[4],
                    target_start=target_interval[0],
                    target_end=target_interval[1],
                    matches=int(fields[9]),
                    block_length=int(fields[10]),
                    mapq=int(fields[11]),
                )
            )
            for field in fields[12:]:
                if field.startswith("pv:Z:"):
                    placement_modes[key].add(field.split(":", 2)[2])

    metrics = {}
    for key in records:
        query = key[1]
        query_union = union_intervals(query_intervals[key])
        target_union = union_intervals(target_intervals[key])
        query_covered_bp = interval_length(query_union)
        target_covered_bp = interval_length(target_union)
        identity = matches[key] / block_lengths[key] if block_lengths[key] else 0.0
        query_coverage = query_covered_bp / query_lengths[query]
        query_span_bp = max(end for start, end in query_intervals[key]) - min(
            start for start, end in query_intervals[key]
        )
        query_span_coverage = query_span_bp / query_lengths[query]
        modes = placement_modes.get(key) or {"dense"}
        if len(modes) != 1:
            raise ValueError(f"Mixed placement modes for {target}/{query}: {sorted(modes)}")
        placement_mode = next(iter(modes))
        placement_support = (
            query_span_coverage if placement_mode == "sparse_syntenic" else query_coverage
        )
        confidence = identity * placement_support * math.log2(1 + target_covered_bp)
        metrics[key] = {
            "identity": identity,
            "query_coverage": query_coverage,
            "query_span_coverage": query_span_coverage,
            "query_covered_bp": query_covered_bp,
            "target_covered_bp": target_covered_bp,
            "confidence": confidence,
            "placement_mode": placement_mode,
        }
    return records, metrics


def merge_projection_blocks(intervals, max_gap):
    blocks = []
    current_intervals = []
    current_end = None
    for alignment in sorted(
        intervals, key=lambda item: (item.target_start, item.target_end)
    ):
        start, end = alignment.target_start, alignment.target_end
        if current_end is None or start <= current_end + max_gap:
            current_intervals.append(alignment)
            current_end = max(current_end or end, end)
        else:
            blocks.append(current_intervals)
            current_intervals = [alignment]
            current_end = end
    if current_intervals:
        blocks.append(current_intervals)
    return blocks


def split_sparse_projection_block(alignments, min_coverage):
    """Split a sparse target block at unsupported gaps until each part is dense."""
    pending = [alignments]
    dense_blocks = []
    while pending:
        candidate = pending.pop()
        block = projection_block_metrics(candidate)
        coverage = block["target_aligned_bp"] / block["target_span_bp"]
        if coverage >= min_coverage or len(candidate) == 1:
            dense_blocks.append(candidate)
            continue

        ordered = sorted(
            candidate, key=lambda item: (item.target_start, item.target_end)
        )
        current_end = ordered[0].target_end
        gaps = []
        for index, alignment in enumerate(ordered[1:], 1):
            if alignment.target_start > current_end:
                gaps.append((alignment.target_start - current_end, index))
            current_end = max(current_end, alignment.target_end)
        if not gaps:
            dense_blocks.append(candidate)
            continue
        _, split_index = max(gaps)
        pending.extend((ordered[:split_index], ordered[split_index:]))
    return sorted(
        dense_blocks,
        key=lambda block: min(alignment.target_start for alignment in block),
    )


def projection_block_metrics(alignments):
    target_intervals = [
        (alignment.target_start, alignment.target_end) for alignment in alignments
    ]
    query_intervals = [
        (alignment.query_start, alignment.query_end) for alignment in alignments
    ]
    target_aligned_bp = interval_length(union_intervals(target_intervals))
    query_aligned_bp = interval_length(union_intervals(query_intervals))
    target_start = min(start for start, _ in target_intervals)
    target_end = max(end for _, end in target_intervals)
    query_start = min(start for start, _ in query_intervals)
    query_end = max(end for _, end in query_intervals)
    block_length = sum(alignment.block_length for alignment in alignments)
    identity = (
        sum(alignment.matches for alignment in alignments) / block_length
        if block_length
        else 0.0
    )
    mean_mapq = (
        sum(alignment.mapq * alignment.block_length for alignment in alignments)
        / block_length
        if block_length
        else 0.0
    )
    strands = {alignment.strand for alignment in alignments}
    strand = next(iter(strands)) if len(strands) == 1 else "mixed"

    ordered = sorted(
        alignments,
        key=lambda item: (
            (item.target_start + item.target_end) / 2,
            item.query_start,
        ),
    )
    oriented_query_midpoints = []
    for alignment in ordered:
        query_midpoint = (alignment.query_start + alignment.query_end) / 2
        if alignment.strand == "-":
            query_midpoint = alignment.query_length - query_midpoint
        oriented_query_midpoints.append(query_midpoint)
    comparisons = len(oriented_query_midpoints) - 1
    collinearity = (
        sum(
            right >= left
            for left, right in zip(
                oriented_query_midpoints, oriented_query_midpoints[1:]
            )
        )
        / comparisons
        if comparisons > 0
        else 1.0
    )
    return {
        "target_start": target_start,
        "target_end": target_end,
        "target_aligned_bp": target_aligned_bp,
        "target_span_bp": target_end - target_start,
        "query_aligned_bp": query_aligned_bp,
        "query_span_bp": query_end - query_start,
        "identity": identity,
        "mean_mapq": mean_mapq,
        "strand": strand,
        "collinearity": collinearity,
    }


def make_projection(
    target,
    unitig,
    block,
    block_index,
    block_count,
    metric,
    dosage,
    contig_type,
    placement_mode,
    args,
    constraint_role="aligned_block",
    projection_class_override=None,
):
    target_fraction = (
        block["target_aligned_bp"] / metric["target_covered_bp"]
        if metric["target_covered_bp"]
        else 0.0
    )
    query_fraction = (
        block["query_aligned_bp"] / metric["query_covered_bp"]
        if metric["query_covered_bp"]
        else 0.0
    )
    projection_coverage = (
        block["target_aligned_bp"] / block["target_span_bp"]
        if block["target_span_bp"]
        else 0.0
    )
    anchor_threshold = max(
        args.min_anchor_block_aligned_bp,
        math.ceil(metric["query_length"] * args.min_anchor_query_fraction),
    )
    projection_class = projection_class_override or (
        "anchor"
        if block["target_aligned_bp"] >= anchor_threshold
        else "supporting"
    )
    if placement_mode == "sparse_syntenic":
        support = metric["query_span_coverage"]
    else:
        support = projection_coverage * math.sqrt(
            max(target_fraction, query_fraction)
        )
    confidence = (
        block["identity"]
        * support
        * math.log2(1 + block["target_aligned_bp"])
    )
    return Projection(
        target=target,
        unitig=unitig,
        start=block["target_start"],
        end=block["target_end"],
        aligned_bp=block["target_aligned_bp"],
        projection_coverage=projection_coverage,
        identity=block["identity"],
        query_coverage=metric["query_coverage"],
        query_span_coverage=metric["query_span_coverage"],
        confidence=confidence,
        dosage=dosage,
        contig_type=contig_type,
        placement_mode=placement_mode,
        block_index=block_index,
        block_count=block_count,
        target_span_bp=block["target_span_bp"],
        query_aligned_bp=block["query_aligned_bp"],
        query_span_bp=block["query_span_bp"],
        target_aligned_fraction=target_fraction,
        query_aligned_fraction=query_fraction,
        chain_identity=metric["identity"],
        mean_mapq=block["mean_mapq"],
        strand=block["strand"],
        collinearity=block["collinearity"],
        projection_class=projection_class,
        query_length=metric["query_length"],
        constraint_role=constraint_role,
    )


def qualifies_for_path_envelope(metric, block, block_count, args):
    """Identify long, chromosome-local paths whose sparse blocks reflect divergence."""
    if metric["placement_mode"] != "segmented":
        return False
    if metric["query_length"] < args.min_path_envelope_query_length:
        return False
    if block_count <= args.max_projection_blocks:
        return False
    if metric["query_coverage"] < args.min_path_envelope_query_coverage:
        return False
    if metric["query_span_coverage"] < args.min_path_envelope_query_span_coverage:
        return False
    if metric["identity"] < args.min_block_identity:
        return False
    if block["mean_mapq"] < args.min_block_mapq:
        return False
    if block["target_span_bp"] <= 0:
        return False
    target_coverage = block["target_aligned_bp"] / block["target_span_bp"]
    if target_coverage < args.min_path_envelope_target_coverage:
        return False
    span_ratio = block["target_span_bp"] / metric["query_length"]
    return span_ratio <= args.max_path_envelope_span_ratio


def build_projections(records, metrics, contig_types, args):
    projections = []
    rejected = []
    fragmented_unitigs = 0

    for (target, unitig), alignments in sorted(records.items()):
        raw_type = contig_types.get(unitig, "unknown")
        if raw_type not in KNOWN_DOSAGE:
            if args.unknown_dosage_policy == "exclude":
                rejected.append((target, unitig, "unknown_dosage_type", raw_type))
                continue
            contig_type = "haplotig"
        else:
            contig_type = raw_type
        dosage = KNOWN_DOSAGE[contig_type]
        metric = metrics[(target, unitig)]
        metric["query_length"] = alignments[0].query_length

        if metric["placement_mode"] == "sparse_syntenic":
            block = projection_block_metrics(alignments)
            if block["target_aligned_bp"] < args.min_projection_aligned_bp:
                rejected.append(
                    (target, unitig, "projection_aligned_bp_below_threshold", "sparse_syntenic")
                )
                continue
            projections.append(
                make_projection(
                    target,
                    unitig,
                    block,
                    1,
                    1,
                    metric,
                    dosage,
                    contig_type,
                    metric["placement_mode"],
                    args,
                )
            )
            continue

        coarse_blocks = merge_projection_blocks(alignments, args.max_projection_gap)
        projection_blocks = [
            block
            for coarse_block in coarse_blocks
            for block in split_sparse_projection_block(
                coarse_block, args.min_projection_coverage
            )
        ]
        if len(projection_blocks) > args.max_projection_blocks:
            fragmented_unitigs += 1

        block_count = len(projection_blocks)
        accepted_block_projections = []
        for block_number, block_alignments in enumerate(projection_blocks, 1):
            block = projection_block_metrics(block_alignments)
            projection_coverage = (
                block["target_aligned_bp"] / block["target_span_bp"]
            )
            detail = f"{block_number}/{block_count}"
            if block["target_aligned_bp"] < args.min_projection_aligned_bp:
                rejected.append(
                    (target, unitig, "projection_aligned_bp_below_threshold", detail)
                )
                continue
            if projection_coverage < args.min_projection_coverage:
                rejected.append(
                    (target, unitig, "projection_coverage_below_threshold", detail)
                )
                continue
            if block["identity"] < args.min_block_identity:
                rejected.append(
                    (target, unitig, "block_identity_below_threshold", detail)
                )
                continue
            if block["mean_mapq"] < args.min_block_mapq:
                rejected.append(
                    (target, unitig, "block_mapq_below_threshold", detail)
                )
                continue
            if (
                block["strand"] == "mixed"
                and metric["placement_mode"] != "segmented"
            ):
                rejected.append((target, unitig, "mixed_block_strands", detail))
                continue
            if (
                block["collinearity"] < args.min_block_collinearity
                and metric["placement_mode"] != "segmented"
            ):
                rejected.append(
                    (target, unitig, "block_collinearity_below_threshold", detail)
                )
                continue
            accepted_block_projections.append(
                make_projection(
                    target,
                    unitig,
                    block,
                    block_number,
                    block_count,
                    metric,
                    dosage,
                    contig_type,
                    metric["placement_mode"],
                    args,
                )
            )
        envelope_block = projection_block_metrics(alignments)
        if accepted_block_projections and qualifies_for_path_envelope(
            metric, envelope_block, block_count, args
        ):
            projections.extend(
                replace(projection, constraint_role="evidence_block")
                for projection in accepted_block_projections
            )
            projections.append(
                make_projection(
                    target,
                    unitig,
                    envelope_block,
                    0,
                    block_count,
                    metric,
                    dosage,
                    contig_type,
                    metric["placement_mode"],
                    args,
                    constraint_role="path_envelope",
                    projection_class_override="path_envelope",
                )
            )
        else:
            projections.extend(accepted_block_projections)
    return projections, rejected, fragmented_unitigs


def top_subsets(candidates, projection_by_unitig, args):
    """Return the two best dosage-valid subsets using confidence and length prior."""
    states = {0: [(0.0, tuple())]}
    for unitig in sorted(candidates):
        projection = projection_by_unitig[unitig]
        additions = defaultdict(list)
        for used, entries in states.items():
            new_used = used + projection.dosage
            if new_used > args.ploidy:
                continue
            for score, subset in entries:
                length_prior = math.sqrt(
                    projection.query_length / args.length_prior_scale
                )
                value = projection.confidence * projection.dosage * length_prior
                additions[new_used].append((score + value, subset + (unitig,)))
        for used, entries in additions.items():
            combined = states.get(used, []) + entries
            unique = {}
            for score, subset in combined:
                unique[subset] = max(score, unique.get(subset, float("-inf")))
            states[used] = sorted(
                ((score, subset) for subset, score in unique.items()), reverse=True
            )[:2]

    ranked = [
        (score, subset)
        for used, entries in states.items()
        if used > 0
        for score, subset in entries
    ]
    return sorted(ranked, reverse=True)[:2]


def resolve_segment(active, gfa_links, args):
    projection_by_unitig = {
        projection.unitig: projection for projection in active
    }
    if len(projection_by_unitig) != len(active):
        raise ValueError("Overlapping projection blocks for the same unitig")
    candidates = tuple(sorted(projection_by_unitig))
    dosage_sum = sum(projection_by_unitig[u].dosage for u in candidates)
    graph_conflicts = tuple(
        (left, right)
        for left, right in combinations(candidates, 2)
        if tuple(sorted((left, right))) in gfa_links
    )
    if graph_conflicts:
        return (
            tuple(),
            tuple(),
            "graph_link_conflict_omitted",
            0.0,
            dosage_sum,
            graph_conflicts,
        )
    if dosage_sum <= args.ploidy:
        return candidates, candidates, "accepted", 1.0, dosage_sum, tuple()

    subsets = top_subsets(candidates, projection_by_unitig, args)
    best_score, best_subset = subsets[0]
    second_score = subsets[1][0] if len(subsets) > 1 else 0.0
    margin = (best_score - second_score) / best_score if best_score else 0.0

    if args.over_capacity_policy == "omit":
        return (
            tuple(),
            tuple(sorted(best_subset)),
            "over_capacity_omitted",
            margin,
            dosage_sum,
            tuple(),
        )
    if args.over_capacity_policy == "confident" and margin < args.min_resolution_margin:
        return (
            tuple(),
            tuple(sorted(best_subset)),
            "over_capacity_ambiguous",
            margin,
            dosage_sum,
            tuple(),
        )
    best_subset = tuple(sorted(best_subset))
    return best_subset, best_subset, "over_capacity_resolved", margin, dosage_sum, tuple()


def sweep_segments(projections, gfa_links, args):
    by_target = defaultdict(list)
    for projection in projections:
        if projection.constraint_role == "evidence_block":
            continue
        by_target[projection.target].append(projection)

    accepted_segments = []
    qc_rows = []
    for target in sorted(by_target):
        events = defaultdict(lambda: {"start": set(), "end": set()})
        for projection in by_target[target]:
            events[projection.start]["start"].add(projection)
            events[projection.end]["end"].add(projection)

        active = set()
        previous = None
        for position in sorted(events):
            if previous is not None and position > previous and active:
                length = position - previous
                if length < args.min_segment_length:
                    selected = tuple()
                    preferred = tuple()
                    status = "segment_too_short"
                    margin = 0.0
                    graph_conflicts = tuple()
                    dosage_sum = sum(
                        projection.dosage for projection in active
                    )
                else:
                    (
                        selected,
                        preferred,
                        status,
                        margin,
                        dosage_sum,
                        graph_conflicts,
                    ) = resolve_segment(active, gfa_links, args)
                    if selected:
                        accepted_segments.append(
                            Segment(target, previous, position, selected)
                        )
                qc_rows.append(
                    {
                        "target": target,
                        "start": previous,
                        "end": position,
                        "length": length,
                        "candidates": ",".join(
                            sorted(projection.unitig for projection in active)
                        ),
                        "candidate_count": len(active),
                        "candidate_dosage": dosage_sum,
                        "selected": ",".join(selected),
                        "preferred": ",".join(preferred),
                        "selected_dosage": sum(
                            next(
                                projection.dosage
                                for projection in active
                                if projection.unitig == unitig
                            )
                            for unitig in selected
                        ),
                        "resolution_margin": margin,
                        "graph_link_conflicts": ",".join(
                            f"{left}|{right}" for left, right in graph_conflicts
                        ),
                        "status": status,
                    }
                )

            # Half-open intervals: ending blocks are inactive after this point;
            # starting blocks are active from this point onward.
            active.difference_update(events[position]["end"])
            active.update(events[position]["start"])
            previous = position

    return accepted_segments, qc_rows


def find_deferred_over_capacity_unitigs(qc_rows, projections, args):
    """Defer short candidates that never occur in a local best dosage subset."""
    projection_by_unitig = {}
    for projection in projections:
        projection_by_unitig.setdefault(projection.unitig, projection)

    candidate_bp = Counter()
    preferred_bp = Counter()
    over_capacity_bp = Counter()
    for row in qc_rows:
        if row["status"] == "accepted" or row["status"].startswith("over_capacity_"):
            candidates = tuple(filter(None, row["candidates"].split(",")))
            preferred = tuple(filter(None, row["preferred"].split(",")))
            for unitig in candidates:
                candidate_bp[unitig] += row["length"]
                if row["status"].startswith("over_capacity_"):
                    over_capacity_bp[unitig] += row["length"]
            for unitig in preferred:
                preferred_bp[unitig] += row["length"]

    deferred = {}
    for unitig, crowded_bp in sorted(over_capacity_bp.items()):
        projection = projection_by_unitig[unitig]
        total_bp = candidate_bp[unitig]
        retained_bp = preferred_bp[unitig]
        preferred_fraction = retained_bp / total_bp if total_bp else 0.0
        if crowded_bp < args.min_deferred_over_capacity_bp:
            continue
        if projection.query_length > args.max_deferred_over_capacity_query_length:
            continue
        if preferred_fraction > args.max_deferred_over_capacity_preferred_fraction:
            continue
        deferred[unitig] = {
            "target": projection.target,
            "unitig": unitig,
            "query_length": projection.query_length,
            "candidate_bp": total_bp,
            "over_capacity_bp": crowded_bp,
            "preferred_bp": retained_bp,
            "preferred_fraction": preferred_fraction,
            "reason": "never_in_preferred_dosage_subset",
        }
    return deferred


def remove_deferred_unitigs(segments, deferred):
    filtered = []
    deferred = set(deferred)
    for segment in segments:
        unitigs = tuple(unitig for unitig in segment.unitigs if unitig not in deferred)
        if unitigs:
            filtered.append(Segment(segment.target, segment.start, segment.end, unitigs))
    return filtered


def subtract_intervals(intervals, covered):
    remaining = []
    covered = union_intervals(covered)
    for start, end in union_intervals(intervals):
        cursor = start
        for covered_start, covered_end in covered:
            if covered_end <= cursor:
                continue
            if covered_start >= end:
                break
            if covered_start > cursor:
                remaining.append((cursor, min(covered_start, end)))
            cursor = max(cursor, covered_end)
            if cursor >= end:
                break
        if cursor < end:
            remaining.append((cursor, end))
    return remaining


def add_long_path_pair_constraints(segments, projections, gfa_links, deferred, args):
    """Preserve long-path versus fragment exclusion independently of ploidy subsets."""
    constraint_projections = [
        projection
        for projection in projections
        if projection.constraint_role != "evidence_block"
        and projection.unitig not in deferred
    ]
    by_target = defaultdict(list)
    for projection in constraint_projections:
        by_target[projection.target].append(projection)

    desired = defaultdict(list)
    skipped_dosage = 0
    skipped_gfa = 0
    for target, target_projections in by_target.items():
        envelopes = [
            projection
            for projection in target_projections
            if projection.constraint_role == "path_envelope"
        ]
        for envelope in envelopes:
            for projection in target_projections:
                if projection.unitig == envelope.unitig:
                    continue
                start = max(envelope.start, projection.start)
                end = min(envelope.end, projection.end)
                if end - start < args.min_long_path_pair_overlap:
                    continue
                pair = tuple(sorted((envelope.unitig, projection.unitig)))
                if pair in gfa_links:
                    skipped_gfa += 1
                    continue
                if envelope.dosage + projection.dosage > args.ploidy:
                    skipped_dosage += 1
                    continue
                desired[(target, pair)].append((start, end))

    covered = defaultdict(list)
    for segment in segments:
        for left, right in combinations(segment.unitigs, 2):
            covered[(segment.target, tuple(sorted((left, right))))].append(
                (segment.start, segment.end)
            )

    supplemental = []
    protected_pairs = set()
    supplemental_bp = 0
    for (target, pair), intervals in sorted(desired.items()):
        protected_pairs.add((target, pair))
        for start, end in subtract_intervals(intervals, covered[(target, pair)]):
            if end <= start:
                continue
            supplemental.append(Segment(target, start, end, pair))
            supplemental_bp += end - start

    metrics = {
        "path_envelopes": sum(
            projection.constraint_role == "path_envelope"
            for projection in constraint_projections
        ),
        "protected_long_path_pairs": len(protected_pairs),
        "supplemental_pair_rows": len(supplemental),
        "supplemental_pair_bp": supplemental_bp,
        "skipped_long_path_pairs_dosage": skipped_dosage,
        "skipped_long_path_pairs_gfa": skipped_gfa,
    }
    return segments + supplemental, metrics


def add_over_capacity_pair_constraints(
    segments, qc_rows, projections, gfa_links, deferred, args
):
    """Recover strong pairwise conflicts hidden by an over-capacity candidate set."""
    projection_by_unitig = {}
    for projection in projections:
        projection_by_unitig.setdefault(projection.unitig, projection)

    desired = defaultdict(list)
    skipped_dosage = 0
    skipped_gfa = 0
    for row in qc_rows:
        if not row["status"].startswith("over_capacity_"):
            continue
        candidates = tuple(
            unitig
            for unitig in filter(None, row["candidates"].split(","))
            if unitig not in deferred
        )
        for left, right in combinations(candidates, 2):
            pair = tuple(sorted((left, right)))
            if pair in gfa_links:
                skipped_gfa += 1
                continue
            if (
                projection_by_unitig[left].dosage
                + projection_by_unitig[right].dosage
                > args.ploidy
            ):
                skipped_dosage += 1
                continue
            desired[(row["target"], pair)].append((row["start"], row["end"]))

    covered = defaultdict(list)
    for segment in segments:
        for left, right in combinations(segment.unitigs, 2):
            covered[(segment.target, tuple(sorted((left, right))))].append(
                (segment.start, segment.end)
            )

    supplemental = []
    recovered_pairs = set()
    supplemental_bp = 0
    below_coverage = 0
    for (target, pair), intervals in sorted(desired.items()):
        intervals = union_intervals(intervals)
        overlap_bp = sum(end - start for start, end in intervals)
        short_length = min(
            projection_by_unitig[pair[0]].query_length,
            projection_by_unitig[pair[1]].query_length,
        )
        short_coverage = overlap_bp / short_length if short_length else 0.0
        if (
            overlap_bp < args.min_over_capacity_pair_overlap
            or short_coverage < args.min_over_capacity_pair_short_coverage
        ):
            below_coverage += 1
            continue
        recovered_pairs.add((target, pair))
        for start, end in subtract_intervals(intervals, covered[(target, pair)]):
            if end <= start:
                continue
            supplemental.append(Segment(target, start, end, pair))
            supplemental_bp += end - start

    metrics = {
        "recovered_pairs": len(recovered_pairs),
        "supplemental_pair_rows": len(supplemental),
        "supplemental_pair_bp": supplemental_bp,
        "below_short_coverage": below_coverage,
        "skipped_pair_dosage": skipped_dosage,
        "skipped_pair_gfa": skipped_gfa,
    }
    return segments + supplemental, metrics


def merge_adjacent_segments(segments):
    merged = []
    for segment in sorted(segments, key=lambda item: (item.target, item.start, item.end)):
        if (
            merged
            and merged[-1].target == segment.target
            and merged[-1].end == segment.start
            and merged[-1].unitigs == segment.unitigs
        ):
            previous = merged[-1]
            merged[-1] = Segment(
                previous.target,
                previous.start,
                segment.end,
                previous.unitigs,
            )
        else:
            merged.append(segment)
    return merged


def write_outputs(
    projections,
    rejected,
    segments,
    qc_rows,
    fragmented_unitigs,
    long_path_metrics,
    over_capacity_pair_metrics,
    deferred_unitigs,
    args,
):
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for segment in segments:
            handle.write(
                f"{segment.target}\t{segment.start}\t{segment.end}\t"
                + "\t".join(segment.unitigs)
                + "\n"
            )

    with args.projections.open("w", newline="") as handle:
        fieldnames = [field.name for field in Projection.__dataclass_fields__.values()]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for projection in projections:
            writer.writerow(projection.__dict__)

    with args.qc.open("w", newline="") as handle:
        fieldnames = list(qc_rows[0]) if qc_rows else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if fieldnames:
            writer.writeheader()
            writer.writerows(qc_rows)

    with args.rejected.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["target", "unitig", "reason", "detail"])
        writer.writerows(rejected)

    deferred_path = args.output.with_name("deferred_overcapacity_unitigs.tsv")
    with deferred_path.open("w", newline="") as handle:
        fieldnames = [
            "target",
            "unitig",
            "query_length",
            "candidate_bp",
            "over_capacity_bp",
            "preferred_bp",
            "preferred_fraction",
            "reason",
        ]
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for unitig in sorted(deferred_unitigs):
            row = dict(deferred_unitigs[unitig])
            row["preferred_fraction"] = f'{row["preferred_fraction"]:.6f}'
            writer.writerow(row)

    unitig_bp = Counter()
    pair_bp = Counter()
    for segment in segments:
        length = segment.end - segment.start
        for unitig in segment.unitigs:
            unitig_bp[(segment.target, unitig)] += length
        for left, right in combinations(segment.unitigs, 2):
            pair_bp[(segment.target, left, right)] += length

    # A path envelope is useful for recovering divergent long paths, but it is
    # inferred across gaps without accepted alignments.  Preserve the overlap
    # between accepted alignment blocks separately so downstream clustering can
    # distinguish direct same-locus evidence from envelope-only evidence.
    direct_intervals = defaultdict(list)
    query_lengths = {}
    for projection in projections:
        query_lengths[projection.unitig] = projection.query_length
        if projection.constraint_role != "path_envelope":
            direct_intervals[(projection.target, projection.unitig)].append(
                (projection.start, projection.end)
            )

    with args.pairs.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "target",
                "unitig1",
                "unitig2",
                "overlap_bp",
                "unitig1_table_bp",
                "unitig2_table_bp",
                "ratio1",
                "ratio2",
                "direct_projection_overlap_bp",
                "direct_query_ratio1",
                "direct_query_ratio2",
                "evidence_class",
            ]
        )
        for (target, left, right), overlap_bp in sorted(pair_bp.items()):
            left_bp = unitig_bp[(target, left)]
            right_bp = unitig_bp[(target, right)]
            direct_overlap_bp = intersection_length(
                direct_intervals[(target, left)],
                direct_intervals[(target, right)],
            )
            left_length = query_lengths[left]
            right_length = query_lengths[right]
            writer.writerow(
                [
                    target,
                    left,
                    right,
                    overlap_bp,
                    left_bp,
                    right_bp,
                    overlap_bp / left_bp,
                    overlap_bp / right_bp,
                    direct_overlap_bp,
                    direct_overlap_bp / left_length if left_length else 0.0,
                    direct_overlap_bp / right_length if right_length else 0.0,
                    (
                        "direct_projection"
                        if direct_overlap_bp >= args.min_segment_length
                        else "boundary_touch"
                        if direct_overlap_bp
                        else "envelope_only"
                    ),
                ]
            )

    summary = {
        "projections": len(projections),
        "rejected_projections": len(rejected),
        "atomic_segments": len(qc_rows),
        "table_rows": len(segments),
        "table_unitigs": len({u for segment in segments for u in segment.unitigs}),
        "projection_mode_counts": dict(Counter(p.placement_mode for p in projections)),
        "projection_class_counts": dict(
            Counter(p.projection_class for p in projections)
        ),
        "projection_role_counts": dict(
            Counter(p.constraint_role for p in projections)
        ),
        "fragmented_unitigs_above_qc_limit": fragmented_unitigs,
        "long_path_constraints": long_path_metrics,
        "over_capacity_pair_constraints": over_capacity_pair_metrics,
        "deferred_over_capacity_unitigs": {
            "count": len(deferred_unitigs),
            "query_bp": sum(
                row["query_length"] for row in deferred_unitigs.values()
            ),
        },
        "status_counts": dict(Counter(row["status"] for row in qc_rows)),
        "parameters": {
            "ploidy": args.ploidy,
            "max_projection_gap": args.max_projection_gap,
            "max_projection_blocks": args.max_projection_blocks,
            "min_projection_aligned_bp": args.min_projection_aligned_bp,
            "min_projection_coverage": args.min_projection_coverage,
            "min_block_identity": args.min_block_identity,
            "min_block_mapq": args.min_block_mapq,
            "min_block_collinearity": args.min_block_collinearity,
            "min_anchor_block_aligned_bp": args.min_anchor_block_aligned_bp,
            "min_anchor_query_fraction": args.min_anchor_query_fraction,
            "min_segment_length": args.min_segment_length,
            "min_path_envelope_query_length": args.min_path_envelope_query_length,
            "min_path_envelope_query_coverage": args.min_path_envelope_query_coverage,
            "min_path_envelope_query_span_coverage": args.min_path_envelope_query_span_coverage,
            "min_path_envelope_target_coverage": args.min_path_envelope_target_coverage,
            "max_path_envelope_span_ratio": args.max_path_envelope_span_ratio,
            "min_long_path_pair_overlap": args.min_long_path_pair_overlap,
            "min_over_capacity_pair_overlap": args.min_over_capacity_pair_overlap,
            "min_over_capacity_pair_short_coverage": args.min_over_capacity_pair_short_coverage,
            "min_deferred_over_capacity_bp": args.min_deferred_over_capacity_bp,
            "max_deferred_over_capacity_query_length": args.max_deferred_over_capacity_query_length,
            "max_deferred_over_capacity_preferred_fraction": args.max_deferred_over_capacity_preferred_fraction,
            "length_prior_scale": args.length_prior_scale,
            "unknown_dosage_policy": args.unknown_dosage_policy,
            "over_capacity_policy": args.over_capacity_policy,
            "min_resolution_margin": args.min_resolution_margin,
            "gfa": str(args.gfa.resolve()) if args.gfa else None,
            "gfa_direct_links": args.gfa_direct_links,
        },
    }
    with args.summary.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Generate a dosage-aware allelic table from interval overlap."
    )
    parser.add_argument("--paf", required=True, type=Path)
    parser.add_argument("--contig-type", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--projections", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path)
    parser.add_argument("--rejected", required=True, type=Path)
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--gfa", type=Path, help="hifiasm GFA used to omit direct-link conflicts")
    parser.add_argument("--ploidy", type=int, default=4)
    parser.add_argument("--max-projection-gap", type=int, default=20000)
    parser.add_argument(
        "--max-projection-blocks",
        type=int,
        default=10,
        help="QC threshold for fragmented unitigs; does not reject the unitig",
    )
    parser.add_argument("--min-projection-aligned-bp", type=int, default=10000)
    parser.add_argument("--min-projection-coverage", type=float, default=0.70)
    parser.add_argument("--min-block-identity", type=float, default=0.90)
    parser.add_argument("--min-block-mapq", type=float, default=20.0)
    parser.add_argument("--min-block-collinearity", type=float, default=0.80)
    parser.add_argument("--min-anchor-block-aligned-bp", type=int, default=100000)
    parser.add_argument("--min-anchor-query-fraction", type=float, default=0.005)
    parser.add_argument("--min-segment-length", type=int, default=10000)
    parser.add_argument("--min-path-envelope-query-length", type=int, default=5000000)
    parser.add_argument("--min-path-envelope-query-coverage", type=float, default=0.25)
    parser.add_argument(
        "--min-path-envelope-query-span-coverage", type=float, default=0.80
    )
    parser.add_argument("--min-path-envelope-target-coverage", type=float, default=0.20)
    parser.add_argument("--max-path-envelope-span-ratio", type=float, default=2.0)
    parser.add_argument("--min-long-path-pair-overlap", type=int, default=10000)
    parser.add_argument("--min-over-capacity-pair-overlap", type=int, default=10000)
    parser.add_argument(
        "--min-over-capacity-pair-short-coverage", type=float, default=0.50
    )
    parser.add_argument("--min-deferred-over-capacity-bp", type=int, default=100000)
    parser.add_argument(
        "--max-deferred-over-capacity-query-length", type=int, default=5000000
    )
    parser.add_argument(
        "--max-deferred-over-capacity-preferred-fraction", type=float, default=0.0
    )
    parser.add_argument("--length-prior-scale", type=float, default=1000000.0)
    parser.add_argument(
        "--unknown-dosage-policy",
        choices=("exclude", "haplotig"),
        default="exclude",
    )
    parser.add_argument(
        "--over-capacity-policy",
        choices=("omit", "confident", "best"),
        default="confident",
    )
    parser.add_argument("--min-resolution-margin", type=float, default=0.10)
    return parser


def main():
    args = build_parser().parse_args()
    if args.min_over_capacity_pair_overlap < 0:
        raise ValueError("--min-over-capacity-pair-overlap must be non-negative")
    if not 0 <= args.min_over_capacity_pair_short_coverage <= 1:
        raise ValueError(
            "--min-over-capacity-pair-short-coverage must be between zero and one"
        )
    if args.min_deferred_over_capacity_bp < 0:
        raise ValueError("--min-deferred-over-capacity-bp must be non-negative")
    if args.max_deferred_over_capacity_query_length < 0:
        raise ValueError(
            "--max-deferred-over-capacity-query-length must be non-negative"
        )
    if not 0 <= args.max_deferred_over_capacity_preferred_fraction <= 1:
        raise ValueError(
            "--max-deferred-over-capacity-preferred-fraction must be between zero and one"
        )
    contig_types = parse_contig_types(args.contig_type)
    gfa_links = parse_gfa_links(args.gfa)
    args.gfa_direct_links = len(gfa_links)
    records, metrics = parse_chain_paf(args.paf)
    projections, rejected, fragmented_unitigs = build_projections(
        records, metrics, contig_types, args
    )
    segments, qc_rows = sweep_segments(projections, gfa_links, args)
    deferred_unitigs = find_deferred_over_capacity_unitigs(
        qc_rows, projections, args
    )
    segments = remove_deferred_unitigs(segments, deferred_unitigs)
    segments = merge_adjacent_segments(segments)
    segments, long_path_metrics = add_long_path_pair_constraints(
        segments, projections, gfa_links, deferred_unitigs, args
    )
    segments = merge_adjacent_segments(segments)
    segments, over_capacity_pair_metrics = add_over_capacity_pair_constraints(
        segments, qc_rows, projections, gfa_links, deferred_unitigs, args
    )
    segments = merge_adjacent_segments(segments)
    write_outputs(
        projections,
        rejected,
        segments,
        qc_rows,
        fragmented_unitigs,
        long_path_metrics,
        over_capacity_pair_metrics,
        deferred_unitigs,
        args,
    )
    print(
        f"Generated {len(segments)} table rows from {len(projections)} projections; "
        f"rejected {len(rejected)} projections"
    )


if __name__ == "__main__":
    main()
