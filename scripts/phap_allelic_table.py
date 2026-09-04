#!/usr/bin/env python3
"""PHap v2 high-confidence allelic unitig table workflow."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def fingerprint(path: Path):
    resolved = path.resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def run_alignment(args, raw_paf: Path):
    if args.paf:
        return args.paf.resolve()
    if raw_paf.exists() and args.reuse_paf:
        return raw_paf
    if raw_paf.exists() and not args.force:
        raise FileExistsError(
            f"{raw_paf} already exists; use --reuse-paf or --force explicitly"
        )

    command = [
        "minimap2",
        "-cx",
        args.minimap_preset,
        "-t",
        str(args.threads),
        str(args.mt2t),
        str(args.p_utg),
    ]
    with raw_paf.open("w") as output:
        subprocess.run(command, stdout=output, check=True)
    return raw_paf


def run_workflow(args):
    script_dir = Path(__file__).resolve().parent
    utils_dir = script_dir.parent / "utils"
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_paf = output_dir / "p_utg_vs_mT2T.paf"
    source_paf = run_alignment(args, raw_paf)
    chain_paf = output_dir / "putg_vs_mT2T.collinear.paf"
    chain_qc = output_dir / "collinear_chain.qc.tsv"
    selection_summary = output_dir / "chromosome_selection.summary.json"

    chain_command = [
        sys.executable,
        str(utils_dir / "find_collinear_chains.py"),
        "--paf",
        str(source_paf),
        "--output",
        str(chain_paf),
        "--qc",
        str(chain_qc),
        "--summary",
        str(selection_summary),
        "--min-alignment-length",
        str(args.min_alignment_length),
        "--min-identity",
        str(args.min_identity),
        "--min-alignment-mapq",
        str(args.min_alignment_mapq),
        "--min-query-length",
        str(args.min_query_length),
        "--min-chain-aligned-bp",
        str(args.min_chain_aligned_bp),
        "--min-query-coverage",
        str(args.min_query_coverage),
        "--min-target-span-coverage",
        str(args.min_target_span_coverage),
        "--min-sparse-query-span-coverage",
        str(args.min_sparse_query_span_coverage),
        "--min-sparse-aligned-bp",
        str(args.min_sparse_aligned_bp),
        "--min-sparse-anchors",
        str(args.min_sparse_anchors),
        "--min-sparse-mean-mapq",
        str(args.min_sparse_mean_mapq),
        "--min-fragmented-aligned-bp",
        str(args.min_fragmented_aligned_bp),
        "--min-fragmented-anchors",
        str(args.min_fragmented_anchors),
        "--min-fragmented-mean-mapq",
        str(args.min_fragmented_mean_mapq),
        "--min-fragmented-reference-margin",
        str(args.min_fragmented_reference_margin),
        "--min-chain-identity",
        str(args.min_chain_identity),
        "--min-reference-margin",
        str(args.min_reference_margin),
        "--max-overlap",
        str(args.max_overlap),
        "--max-gap",
        str(args.max_gap),
    ]
    subprocess.run(chain_command, check=True)

    table_path = output_dir / "corrected_allelic_table.txt"
    table_alias = output_dir / "allelic.ctg.table.v2"
    projections = output_dir / "unitig_projections.tsv"
    table_qc = output_dir / "allelic_table.qc.tsv"
    rejected = output_dir / "rejected_projections.tsv"
    pairs = output_dir / "allelic_pairs.tsv"
    summary = output_dir / "allelic_table.summary.json"

    table_command = [
        sys.executable,
        str(utils_dir / "allelic_table_generate_v2.py"),
        "--paf",
        str(chain_paf),
        "--output",
        str(table_path),
        "--projections",
        str(projections),
        "--qc",
        str(table_qc),
        "--rejected",
        str(rejected),
        "--pairs",
        str(pairs),
        "--summary",
        str(summary),
        "--ploidy",
        str(args.ploidy),
        "--max-projection-gap",
        str(args.max_projection_gap),
        "--max-projection-blocks",
        str(args.max_projection_blocks),
        "--min-projection-aligned-bp",
        str(args.min_projection_aligned_bp),
        "--min-projection-coverage",
        str(args.min_projection_coverage),
        "--min-block-identity",
        str(args.min_block_identity),
        "--min-block-mapq",
        str(args.min_block_mapq),
        "--min-block-collinearity",
        str(args.min_block_collinearity),
        "--min-anchor-block-aligned-bp",
        str(args.min_anchor_block_aligned_bp),
        "--min-anchor-query-fraction",
        str(args.min_anchor_query_fraction),
        "--min-segment-length",
        str(args.min_segment_length),
        "--min-path-envelope-query-length",
        str(args.min_path_envelope_query_length),
        "--min-path-envelope-query-coverage",
        str(args.min_path_envelope_query_coverage),
        "--min-path-envelope-query-span-coverage",
        str(args.min_path_envelope_query_span_coverage),
        "--min-path-envelope-target-coverage",
        str(args.min_path_envelope_target_coverage),
        "--max-path-envelope-span-ratio",
        str(args.max_path_envelope_span_ratio),
        "--min-long-path-pair-overlap",
        str(args.min_long_path_pair_overlap),
        "--min-over-capacity-pair-overlap",
        str(args.min_over_capacity_pair_overlap),
        "--min-over-capacity-pair-short-coverage",
        str(args.min_over_capacity_pair_short_coverage),
        "--min-deferred-over-capacity-bp",
        str(args.min_deferred_over_capacity_bp),
        "--max-deferred-over-capacity-query-length",
        str(args.max_deferred_over_capacity_query_length),
        "--max-deferred-over-capacity-preferred-fraction",
        str(args.max_deferred_over_capacity_preferred_fraction),
        "--length-prior-scale",
        str(args.length_prior_scale),
        "--unknown-dosage-policy",
        args.unknown_dosage_policy,
        "--over-capacity-policy",
        args.over_capacity_policy,
        "--min-resolution-margin",
        str(args.min_resolution_margin),
    ]
    if args.no_collapse:
        table_command.append("--no-collapse")
    else:
        table_command.extend(["--contig-type", str(args.contig_type)])
    if args.gfa:
        table_command.extend(["--gfa", str(args.gfa)])
    subprocess.run(table_command, check=True)
    shutil.copyfile(table_path, table_alias)

    gfa_validation = output_dir / "allelic_pairs.gfa_validation.tsv"
    gfa_summary = output_dir / "allelic_pairs.gfa_validation.summary.json"
    if args.gfa:
        subprocess.run(
            [
                sys.executable,
                str(utils_dir / "validate_allelic_table_gfa.py"),
                "--gfa",
                str(args.gfa),
                "--pairs",
                str(pairs),
                "--output",
                str(gfa_validation),
                "--summary",
                str(gfa_summary),
            ],
            check=True,
        )

    manifest = {
        "phap_version": "2.0.0",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "p_utg": fingerprint(args.p_utg),
            "mt2t": fingerprint(args.mt2t),
            "contig_type": (
                fingerprint(args.contig_type) if args.contig_type else None
            ),
            "paf": fingerprint(source_paf),
        },
        "outputs": {
            "collinear_paf": str(chain_paf),
            "chain_qc": str(chain_qc),
            "selection_summary": str(selection_summary),
            "allelic_table": str(table_path),
            "table_qc": str(table_qc),
            "projections": str(projections),
            "pairs": str(pairs),
            "summary": str(summary),
            "deferred_overcapacity": str(
                output_dir / "deferred_overcapacity_unitigs.tsv"
            ),
        },
        "parameters": vars(args) | {"output_dir": str(args.output_dir)},
    }
    if args.gfa:
        manifest["inputs"]["gfa"] = fingerprint(args.gfa)
        manifest["outputs"]["gfa_validation"] = str(gfa_validation)
        manifest["outputs"]["gfa_validation_summary"] = str(gfa_summary)
    for key, value in list(manifest["parameters"].items()):
        if isinstance(value, Path):
            manifest["parameters"][key] = str(value.resolve())
    with (output_dir / "run_manifest.json").open("w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"PHap v2 allelic table: {table_path}")
    print(f"QC summary: {summary}")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Generate a high-confidence PHap v2 allelic unitig table."
    )
    parser.add_argument("--p_utg", required=True, type=Path)
    parser.add_argument("--mT2T", dest="mt2t", required=True, type=Path)
    parser.add_argument(
        "--contig_type", type=Path,
        help="Contig dosage table; omitted with --no-collapse",
    )
    parser.add_argument(
        "--no-collapse", "--no_collapse", dest="no_collapse", action="store_true",
        help="Treat every unitig as single-copy dosage one",
    )
    parser.add_argument("--paf", type=Path, help="Reuse an explicitly supplied raw PAF")
    parser.add_argument("--gfa", type=Path, help="hifiasm GFA for graph-aware validation")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("02.cluster.v2/01.putg_vs_mT2T"),
    )
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--minimap-preset", default="asm5")
    parser.add_argument("--reuse-paf", action="store_true")
    parser.add_argument("--force", action="store_true")

    chain = parser.add_argument_group("chromosome-local alignment filters")
    chain.add_argument("--min-alignment-length", type=int, default=1000)
    chain.add_argument("--min-identity", type=float, default=0.90)
    chain.add_argument("--min-alignment-mapq", type=int, default=20)
    chain.add_argument("--min-query-length", type=int, default=20000)
    chain.add_argument("--min-chain-aligned-bp", type=int, default=20000)
    chain.add_argument("--min-query-coverage", type=float, default=0.30)
    chain.add_argument("--min-target-span-coverage", type=float, default=0.30)
    chain.add_argument("--min-sparse-query-span-coverage", type=float, default=0.80)
    chain.add_argument("--min-sparse-aligned-bp", type=int, default=1000000)
    chain.add_argument("--min-sparse-anchors", type=int, default=20)
    chain.add_argument("--min-sparse-mean-mapq", type=float, default=20.0)
    chain.add_argument("--min-fragmented-aligned-bp", type=int, default=500000)
    chain.add_argument("--min-fragmented-anchors", type=int, default=20)
    chain.add_argument("--min-fragmented-mean-mapq", type=float, default=20.0)
    chain.add_argument("--min-fragmented-reference-margin", type=float, default=0.20)
    chain.add_argument("--min-chain-identity", type=float, default=0.90)
    chain.add_argument("--min-reference-margin", type=float, default=0.05)
    chain.add_argument("--max-overlap", type=int, default=10000)
    chain.add_argument("--max-gap", type=int, default=3000000)

    table = parser.add_argument_group("allelic table filters")
    table.add_argument("--ploidy", type=int, default=4)
    table.add_argument("--max-projection-gap", type=int, default=20000)
    table.add_argument(
        "--max-projection-blocks",
        type=int,
        default=10,
        help="QC threshold for fragmented unitigs; does not reject the unitig",
    )
    table.add_argument("--min-projection-aligned-bp", type=int, default=10000)
    table.add_argument("--min-projection-coverage", type=float, default=0.70)
    table.add_argument("--min-block-identity", type=float, default=0.90)
    table.add_argument("--min-block-mapq", type=float, default=20.0)
    table.add_argument("--min-block-collinearity", type=float, default=0.80)
    table.add_argument("--min-anchor-block-aligned-bp", type=int, default=100000)
    table.add_argument("--min-anchor-query-fraction", type=float, default=0.005)
    table.add_argument("--min-segment-length", type=int, default=10000)
    table.add_argument("--min-path-envelope-query-length", type=int, default=5000000)
    table.add_argument("--min-path-envelope-query-coverage", type=float, default=0.25)
    table.add_argument(
        "--min-path-envelope-query-span-coverage", type=float, default=0.80
    )
    table.add_argument("--min-path-envelope-target-coverage", type=float, default=0.20)
    table.add_argument("--max-path-envelope-span-ratio", type=float, default=2.0)
    table.add_argument("--min-long-path-pair-overlap", type=int, default=10000)
    table.add_argument("--min-over-capacity-pair-overlap", type=int, default=10000)
    table.add_argument(
        "--min-over-capacity-pair-short-coverage", type=float, default=0.50
    )
    table.add_argument("--min-deferred-over-capacity-bp", type=int, default=100000)
    table.add_argument(
        "--max-deferred-over-capacity-query-length", type=int, default=5000000
    )
    table.add_argument(
        "--max-deferred-over-capacity-preferred-fraction", type=float, default=0.0
    )
    table.add_argument("--length-prior-scale", type=float, default=1000000.0)
    table.add_argument(
        "--unknown-dosage-policy",
        choices=("exclude", "haplotig"),
        default="exclude",
    )
    table.add_argument(
        "--over-capacity-policy",
        choices=("omit", "confident", "best"),
        default="confident",
    )
    table.add_argument("--min-resolution-margin", type=float, default=0.10)
    return parser


def main():
    args = build_parser().parse_args()
    if args.no_collapse and args.contig_type is not None:
        raise ValueError("--no-collapse and --contig_type are mutually exclusive")
    if not args.no_collapse and args.contig_type is None:
        raise ValueError("--contig_type is required unless --no-collapse is used")
    run_workflow(args)


if __name__ == "__main__":
    main()
