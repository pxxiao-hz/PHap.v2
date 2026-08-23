#!/usr/bin/env python3
"""Cross-check predicted allelic pairs against a hifiasm GFA graph."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict, deque
from pathlib import Path


def parse_pairs(path: Path):
    rows = []
    wanted = set()
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)
            wanted.update((row["unitig1"], row["unitig2"]))
    return rows, wanted


def parse_overlap(cigar: str) -> int:
    match = re.fullmatch(r"(\d+)M", cigar)
    return int(match.group(1)) if match else 0


def parse_gfa(path: Path, wanted):
    lengths = {}
    adjacency = defaultdict(set)
    link_overlap = {}
    reads = defaultdict(set)

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            if fields[0] == "S":
                length_tags = [tag for tag in fields[3:] if tag.startswith("LN:i:")]
                if length_tags:
                    lengths[fields[1]] = int(length_tags[0][5:])
            elif fields[0] == "L":
                if len(fields) < 6:
                    raise ValueError(f"Malformed GFA L record at line {line_number}")
                left, right = fields[1], fields[3]
                key = tuple(sorted((left, right)))
                adjacency[left].add(right)
                adjacency[right].add(left)
                link_overlap[key] = max(
                    link_overlap.get(key, 0), parse_overlap(fields[5])
                )
            elif fields[0] == "A" and fields[1] in wanted:
                if len(fields) < 5:
                    raise ValueError(f"Malformed GFA A record at line {line_number}")
                reads[fields[1]].add(fields[4])

    return lengths, adjacency, link_overlap, reads


def capped_distance(adjacency, source, target, limit):
    if source == target:
        return 0
    visited = {source}
    queue = deque([(source, 0)])
    while queue:
        node, distance = queue.popleft()
        if distance == limit:
            continue
        for neighbor in adjacency.get(node, ()):
            if neighbor == target:
                return distance + 1
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, distance + 1))
    return None


def classify(direct_link, common_neighbor_count, shared_read_fraction):
    if direct_link:
        return "conflict_direct_link"
    if common_neighbor_count >= 2:
        return "support_bubble_like"
    if shared_read_fraction >= 0.05:
        return "review_shared_reads"
    return "graph_neutral"


def validate(rows, lengths, adjacency, link_overlap, reads, max_distance):
    output_rows = []
    for row in rows:
        left, right = row["unitig1"], row["unitig2"]
        key = tuple(sorted((left, right)))
        common_neighbors = sorted(adjacency[left] & adjacency[right])
        shared_reads = len(reads[left] & reads[right])
        minimum_reads = min(len(reads[left]), len(reads[right]))
        shared_fraction = shared_reads / minimum_reads if minimum_reads else 0.0
        graph_distance = capped_distance(adjacency, left, right, max_distance)
        direct_link = key in link_overlap
        gfa_overlap = link_overlap.get(key, 0)
        table_overlap = int(row["overlap_bp"])
        output_rows.append(
            {
                **row,
                "unitig1_gfa_length": lengths.get(left, ""),
                "unitig2_gfa_length": lengths.get(right, ""),
                "unitig1_degree": len(adjacency[left]),
                "unitig2_degree": len(adjacency[right]),
                "direct_gfa_link": int(direct_link),
                "gfa_link_overlap_bp": gfa_overlap,
                "table_to_gfa_overlap_ratio": (
                    table_overlap / gfa_overlap if gfa_overlap else ""
                ),
                "common_neighbor_count": len(common_neighbors),
                "common_neighbors": ",".join(common_neighbors),
                "graph_distance": graph_distance if graph_distance is not None else f">{max_distance}",
                "unitig1_gfa_reads": len(reads[left]),
                "unitig2_gfa_reads": len(reads[right]),
                "shared_gfa_reads": shared_reads,
                "shared_read_fraction_min": shared_fraction,
                "gfa_evidence": classify(
                    direct_link, len(common_neighbors), shared_fraction
                ),
            }
        )
    return output_rows


def write_outputs(rows, output: Path, summary_path: Path, gfa: Path, max_distance):
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else []
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if fieldnames:
            writer.writeheader()
            writer.writerows(rows)

    evidence_counts = Counter(row["gfa_evidence"] for row in rows)
    summary = {
        "gfa": str(gfa.resolve()),
        "pairs": len(rows),
        "gfa_evidence_counts": dict(evidence_counts),
        "direct_link_pairs": sum(int(row["direct_gfa_link"]) for row in rows),
        "pairs_with_shared_reads": sum(int(row["shared_gfa_reads"]) > 0 for row in rows),
        "pairs_with_two_or_more_common_neighbors": sum(
            int(row["common_neighbor_count"]) >= 2 for row in rows
        ),
        "max_graph_distance": max_distance,
    }
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Validate allelic pairs against hifiasm GFA links and reads."
    )
    parser.add_argument("--gfa", required=True, type=Path)
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--max-distance", type=int, default=4)
    return parser


def main():
    args = build_parser().parse_args()
    pair_rows, wanted = parse_pairs(args.pairs)
    graph = parse_gfa(args.gfa, wanted)
    rows = validate(pair_rows, *graph, args.max_distance)
    write_outputs(rows, args.output, args.summary, args.gfa, args.max_distance)
    counts = Counter(row["gfa_evidence"] for row in rows)
    print(f"Validated {len(rows)} pairs: {dict(counts)}")


if __name__ == "__main__":
    main()
