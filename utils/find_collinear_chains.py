#!/usr/bin/env python3
"""Select chromosome-specific local PAF evidence for each query unitig.

Autopolyploid p_utg sequences can contain collapsed sequence, inversions, and
haplotype-specific rearrangements. Chromosome assignment therefore aggregates
reliable local alignments on each target without requiring one strand or one
globally monotonic chain. All retained blocks on the winning chromosome are
written for interval-level projection by the allelic-table stage.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Alignment:
    fields: tuple[str, ...]
    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target: str
    target_length: int
    target_start: int
    target_end: int
    matches: int
    block_length: int
    mapq: int

    @property
    def identity(self) -> float:
        return self.matches / self.block_length if self.block_length else 0.0


@dataclass
class Placement:
    query: str
    target: str
    alignments: list[Alignment]
    score: float
    query_covered_bp: int
    query_span_bp: int
    target_covered_bp: int
    target_span_bp: int
    matches: int
    block_length: int
    mean_mapq: float
    strand: str
    strand_switches: int
    placement_mode: str = "segmented"

    @property
    def identity(self) -> float:
        return self.matches / self.block_length if self.block_length else 0.0

    @property
    def query_coverage(self) -> float:
        query_length = self.alignments[0].query_length
        return self.query_covered_bp / query_length if query_length else 0.0

    @property
    def query_span_coverage(self) -> float:
        query_length = self.alignments[0].query_length
        return self.query_span_bp / query_length if query_length else 0.0

    @property
    def target_span_coverage(self) -> float:
        return self.target_covered_bp / self.target_span_bp if self.target_span_bp else 0.0


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_paf(
    path: Path,
    min_alignment_length: int,
    min_identity: float,
    min_alignment_mapq: int,
):
    grouped: dict[str, dict[str, list[Alignment]]] = defaultdict(
        lambda: defaultdict(list)
    )
    query_lengths: dict[str, int] = {}
    raw_counts = defaultdict(int)
    filtered_counts = defaultdict(int)

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF line {line_number}: fewer than 12 fields")
            tags = parse_tags(fields)
            query = fields[0]
            query_length = int(fields[1])
            previous_length = query_lengths.setdefault(query, query_length)
            if previous_length != query_length:
                raise ValueError(f"Inconsistent query length for {query}")
            raw_counts[query] += 1

            # Secondary hits are repeat alternatives and must not vote for a
            # chromosome. Supplementary primary blocks remain valid evidence.
            if tags.get("tp", "P") != "P":
                filtered_counts[query] += 1
                continue
            alignment = Alignment(
                fields=tuple(fields),
                query=query,
                query_length=query_length,
                query_start=int(fields[2]),
                query_end=int(fields[3]),
                strand=fields[4],
                target=fields[5],
                target_length=int(fields[6]),
                target_start=int(fields[7]),
                target_end=int(fields[8]),
                matches=int(fields[9]),
                block_length=int(fields[10]),
                mapq=int(fields[11]),
            )
            if (
                alignment.block_length < min_alignment_length
                or alignment.identity < min_identity
                or alignment.mapq < min_alignment_mapq
            ):
                filtered_counts[query] += 1
                continue
            grouped[query][alignment.target].append(alignment)

    return grouped, query_lengths, raw_counts, filtered_counts


def union_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def weighted_query_evidence(alignments: list[Alignment]) -> float:
    """Integrate the strongest alignment quality at each query position."""
    events = defaultdict(lambda: {"start": [], "end": []})
    active = [False] * len(alignments)
    for index, alignment in enumerate(alignments):
        mapq_weight = 0.5 + 0.5 * min(alignment.mapq, 60) / 60
        weight = alignment.identity * mapq_weight
        events[alignment.query_start]["start"].append((index, weight))
        events[alignment.query_end]["end"].append(index)

    heap = []
    score = 0.0
    previous = None
    for position in sorted(events):
        while heap and not active[heap[0][1]]:
            heapq.heappop(heap)
        if previous is not None and position > previous and heap:
            score += (position - previous) * -heap[0][0]
        for index in events[position]["end"]:
            active[index] = False
        for index, weight in events[position]["start"]:
            active[index] = True
            heapq.heappush(heap, (-weight, index))
        previous = position
    return score


def summarize_placement(alignments: list[Alignment]) -> Placement:
    query_union = union_intervals(
        (item.query_start, item.query_end) for item in alignments
    )
    target_union = union_intervals(
        (item.target_start, item.target_end) for item in alignments
    )
    strands = {item.strand for item in alignments}
    strand = next(iter(strands)) if len(strands) == 1 else "mixed"
    ordered = sorted(
        alignments,
        key=lambda item: (item.query_start, item.query_end, item.target_start),
    )
    strand_switches = sum(
        left.strand != right.strand for left, right in zip(ordered, ordered[1:])
    )
    block_length = sum(item.block_length for item in alignments)
    return Placement(
        query=alignments[0].query,
        target=alignments[0].target,
        alignments=alignments,
        score=weighted_query_evidence(alignments),
        query_covered_bp=sum(end - start for start, end in query_union),
        query_span_bp=query_union[-1][1] - query_union[0][0],
        target_covered_bp=sum(end - start for start, end in target_union),
        target_span_bp=target_union[-1][1] - target_union[0][0],
        matches=sum(item.matches for item in alignments),
        block_length=block_length,
        mean_mapq=(
            sum(item.mapq * item.block_length for item in alignments) / block_length
            if block_length
            else 0.0
        ),
        strand=strand,
        strand_switches=strand_switches,
    )


def select_placements(grouped, query_lengths, raw_counts, filtered_counts, args):
    accepted: list[Placement] = []
    qc_rows = []

    for query in sorted(query_lengths):
        query_length = query_lengths[query]
        candidates = [
            summarize_placement(alignments)
            for alignments in grouped.get(query, {}).values()
            if alignments
        ]
        candidates.sort(
            key=lambda item: (item.score, item.query_covered_bp, item.identity, item.target),
            reverse=True,
        )
        best = candidates[0] if candidates else None
        second = candidates[1] if len(candidates) > 1 else None
        margin = 1.0
        if best and second and best.score > 0:
            margin = (best.score - second.score) / best.score

        status = "accepted"
        reason = ""
        if query_length < args.min_query_length:
            status, reason = "rejected", "query_too_short"
        elif best is None:
            status, reason = "rejected", "no_primary_alignment_after_filters"
        elif best.query_covered_bp < args.min_chain_aligned_bp:
            status, reason = "rejected", "target_evidence_bp_below_threshold"
        elif best.identity < args.min_chain_identity:
            status, reason = "rejected", "target_evidence_identity_below_threshold"
        elif second is not None and margin < args.min_reference_margin:
            status, reason = "rejected", "ambiguous_reference_chromosome"

        if status == "accepted":
            accepted.append(best)

        qc_rows.append(
            {
                "query": query,
                "query_length": query_length,
                "raw_alignments": raw_counts.get(query, 0),
                "filtered_alignments": filtered_counts.get(query, 0),
                "candidate_targets": len(candidates),
                "best_target": best.target if best else "",
                "best_strand": best.strand if best else "",
                "selected_alignments": len(best.alignments) if best else 0,
                "chain_alignments": len(best.alignments) if best else 0,
                "query_covered_bp": best.query_covered_bp if best else 0,
                "query_coverage": best.query_coverage if best else 0,
                "query_span_bp": best.query_span_bp if best else 0,
                "query_span_coverage": best.query_span_coverage if best else 0,
                "target_covered_bp": best.target_covered_bp if best else 0,
                "target_span_bp": best.target_span_bp if best else 0,
                "target_span_coverage": best.target_span_coverage if best else 0,
                "chain_identity": best.identity if best else 0,
                "mean_mapq": best.mean_mapq if best else 0,
                "chain_score": best.score if best else 0,
                "second_score": second.score if second else 0,
                "reference_margin": margin,
                "strand_switches": best.strand_switches if best else 0,
                "placement_mode": best.placement_mode if status == "accepted" else "",
                "status": status,
                "reason": reason,
            }
        )
    return accepted, qc_rows


def write_outputs(placements, qc_rows, output_paf: Path, qc_path: Path):
    output_paf.parent.mkdir(parents=True, exist_ok=True)
    with output_paf.open("w") as handle:
        for placement in sorted(placements, key=lambda item: (item.query, item.target)):
            for alignment in sorted(
                placement.alignments,
                key=lambda item: (
                    item.query_start,
                    item.query_end,
                    item.target_start,
                    item.target_end,
                    item.strand,
                ),
            ):
                fields = [field for field in alignment.fields if not field.startswith("pv:")]
                fields.append(f"pv:Z:{placement.placement_mode}")
                handle.write("\t".join(fields) + "\n")

    fieldnames = list(qc_rows[0]) if qc_rows else []
    with qc_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if fieldnames:
            writer.writeheader()
            writer.writerows(qc_rows)


def write_summary(path: Path, placements, qc_rows, args):
    status_counts = Counter(row["status"] for row in qc_rows)
    reason_counts = Counter(
        row["reason"] for row in qc_rows if row["reason"]
    )
    summary = {
        "queries": len(qc_rows),
        "accepted_queries": len(placements),
        "accepted_query_bp": sum(
            placement.alignments[0].query_length for placement in placements
        ),
        "selected_paf_records": sum(
            len(placement.alignments) for placement in placements
        ),
        "mixed_strand_queries": sum(
            placement.strand == "mixed" for placement in placements
        ),
        "queries_with_strand_switches": sum(
            placement.strand_switches > 0 for placement in placements
        ),
        "status_counts": dict(sorted(status_counts.items())),
        "rejection_reason_counts": dict(sorted(reason_counts.items())),
        "parameters": {
            "min_alignment_length": args.min_alignment_length,
            "min_identity": args.min_identity,
            "min_alignment_mapq": args.min_alignment_mapq,
            "min_query_length": args.min_query_length,
            "min_target_evidence_bp": args.min_chain_aligned_bp,
            "min_target_evidence_identity": args.min_chain_identity,
            "min_reference_margin": args.min_reference_margin,
            "global_collinearity_required": False,
            "single_strand_required": False,
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Select chromosome-specific, potentially rearranged local PAF evidence."
    )
    parser.add_argument("--paf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--qc", required=True, type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--min-alignment-length", type=int, default=1000)
    parser.add_argument("--min-identity", type=float, default=0.90)
    parser.add_argument("--min-alignment-mapq", type=int, default=20)
    parser.add_argument("--min-query-length", type=int, default=20000)
    parser.add_argument("--min-chain-aligned-bp", type=int, default=20000)
    parser.add_argument("--min-chain-identity", type=float, default=0.90)
    parser.add_argument("--min-reference-margin", type=float, default=0.05)

    # Kept for command-line compatibility. Global collinearity is not assumed.
    parser.add_argument("--min-query-coverage", type=float, default=0.30)
    parser.add_argument("--min-target-span-coverage", type=float, default=0.30)
    parser.add_argument("--min-sparse-query-span-coverage", type=float, default=0.80)
    parser.add_argument("--min-sparse-aligned-bp", type=int, default=1000000)
    parser.add_argument("--min-sparse-anchors", type=int, default=20)
    parser.add_argument("--min-sparse-mean-mapq", type=float, default=20.0)
    parser.add_argument("--min-fragmented-aligned-bp", type=int, default=500000)
    parser.add_argument("--min-fragmented-anchors", type=int, default=20)
    parser.add_argument("--min-fragmented-mean-mapq", type=float, default=20.0)
    parser.add_argument("--min-fragmented-reference-margin", type=float, default=0.20)
    parser.add_argument("--max-overlap", type=int, default=10000)
    parser.add_argument("--max-gap", type=int, default=3000000)
    parser.add_argument("--overlap-penalty", type=float, default=1.0)
    parser.add_argument("--gap-discordance-penalty", type=float, default=0.01)
    return parser


def main():
    args = build_parser().parse_args()
    grouped, query_lengths, raw_counts, filtered_counts = parse_paf(
        args.paf,
        args.min_alignment_length,
        args.min_identity,
        args.min_alignment_mapq,
    )
    placements, qc_rows = select_placements(
        grouped, query_lengths, raw_counts, filtered_counts, args
    )
    write_outputs(placements, qc_rows, args.output, args.qc)
    if args.summary:
        write_summary(args.summary, placements, qc_rows, args)
    print(f"Accepted {len(placements)} of {len(query_lengths)} query unitigs")


if __name__ == "__main__":
    main()
