#!/usr/bin/env python3
"""Materialize conflict-aware mT2T overlap paths and complete audit tables."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from typing import Mapping, Optional

from phap_core import __version__
from phap_core.atomic_io import atomic_write_lines
from phap_core.fasta import FastaRecord, fasta_lines, read_fasta
from phap_core.mt2t_overlap import Mt2tAssemblyResult, Mt2tOverlapParameters, assemble_mt2t
from phap_core.paf import read_paf


def run_mt2t_overlap_workflow(
    *,
    fasta_path: str,
    paf_path: str,
    output_directory: str,
    parameters: Mt2tOverlapParameters,
    tool_versions: Optional[Mapping[str, str]] = None,
) -> Mt2tAssemblyResult:
    """Run the pure graph algorithm and atomically write its products."""

    output_root = Path(output_directory)
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = output_root / "mt2t_manifest.tsv"
    manifest_path.unlink(missing_ok=True)
    fasta_records = read_fasta(fasta_path)
    result = assemble_mt2t(fasta_records, read_paf(paf_path), parameters)
    records_by_id = {record.identifier: record for record in fasta_records}
    _validate_routing(result, records_by_id)

    atomic_write_lines(
        output_root / "mT2T.fa",
        (
            line
            for output in result.outputs
            for line in (f">{output.header}", output.sequence)
        ),
    )
    retained = tuple(
        record for record in fasta_records if record.identifier not in result.contained_ids
    )
    removed = tuple(
        record for record in fasta_records if record.identifier in result.contained_ids
    )
    atomic_write_lines(output_root / "retained.fa", fasta_lines(retained))
    atomic_write_lines(output_root / "removed.fa", fasta_lines(removed))

    atomic_write_lines(
        output_root / "mt2t_alignment_audit.tsv",
        (
            "source\tline\tquery_ID\ttarget_ID\tinitial_status\tinitial_reason\t"
            "final_status\tfinal_reason",
            *(
                "\t".join(
                    (
                        row.source,
                        str(row.line_number),
                        row.query_id,
                        row.target_id,
                        row.initial_status,
                        row.initial_reason,
                        row.final_status,
                        row.final_reason,
                    )
                )
                for row in result.alignment_audit
            ),
        ),
    )
    atomic_write_lines(
        output_root / "mt2t_chain_audit.tsv",
        (
            "query_ID\ttarget_ID\tstrand\trecord_count\tidentity\tquery_coverage\t"
            "target_coverage\tquery_union_bases\ttarget_union_bases\t"
            "query_span_bases\ttarget_span_bases\tscore\tcollinear\t"
            "classification\treason",
            *(
                "\t".join(
                    (
                        row.query_id,
                        row.target_id,
                        row.strand,
                        str(row.record_count),
                        f"{row.identity:.12g}",
                        f"{row.query_coverage:.12g}",
                        f"{row.target_coverage:.12g}",
                        str(row.query_union_bases),
                        str(row.target_union_bases),
                        str(row.query_span_bases),
                        str(row.target_span_bases),
                        f"{row.score:.12g}",
                        str(row.collinear).lower(),
                        row.classification,
                        row.reason,
                    )
                )
                for row in result.chain_audit
            ),
        ),
    )
    atomic_write_lines(
        output_root / "mt2t_overlap_edges.tsv",
        (
            "left_ID\tleft_orientation\tleft_endpoint\tright_ID\tright_orientation\t"
            "right_endpoint\tleft_start\tleft_end\tright_start\tright_end\tidentity\t"
            "query_coverage\ttarget_coverage\tscore\tstatus\treason\tsource_records",
            *(
                "\t".join(
                    (
                        row.edge.left.contig_id,
                        row.edge.left.orientation,
                        row.edge.left_endpoint[1],
                        row.edge.right.contig_id,
                        row.edge.right.orientation,
                        row.edge.right_endpoint[1],
                        str(row.edge.left_interval[0]),
                        str(row.edge.left_interval[1]),
                        str(row.edge.right_interval[0]),
                        str(row.edge.right_interval[1]),
                        f"{row.edge.identity:.12g}",
                        f"{row.edge.query_coverage:.12g}",
                        f"{row.edge.target_coverage:.12g}",
                        f"{row.edge.score:.12g}",
                        row.status,
                        row.reason,
                        ",".join(
                            f"{source}:{line}"
                            for source, line in row.edge.source_records
                        ),
                    )
                )
                for row in result.edge_audit
            ),
        ),
    )
    atomic_write_lines(
        output_root / "mt2t_containment_audit.tsv",
        (
            "source_ID\thost_ID\tstrand\tidentity\tchild_coverage\t"
            "host_start\thost_end\tsource_records\tscore\tstatus\treason",
            *(
                "\t".join(
                    (
                        row.child_id,
                        row.host_id,
                        row.strand,
                        f"{row.identity:.12g}",
                        f"{row.child_coverage:.12g}",
                        str(row.host_start),
                        str(row.host_end),
                        ",".join(
                            f"{source}:{line}"
                            for source, line in row.source_records
                        ),
                        f"{row.score:.12g}",
                        row.status,
                        row.reason,
                    )
                )
                for row in result.containment_audit
            ),
        ),
    )
    atomic_write_lines(
        output_root / "mt2t_join_audit.tsv",
        (
            "destination_ID\tstep\tleft_ID\tleft_orientation\tright_ID\t"
            "right_orientation\tleft_start\tleft_end\tright_start\tright_end\t"
            "right_trim_bases\tright_retained_bases\tscore",
            *(
                "\t".join(
                    (
                        row.destination_id,
                        str(row.step),
                        row.left_id,
                        row.left_orientation,
                        row.right_id,
                        row.right_orientation,
                        str(row.left_start),
                        str(row.left_end),
                        str(row.right_start),
                        str(row.right_end),
                        str(row.right_trim_bases),
                        str(row.right_retained_bases),
                        f"{row.score:.12g}",
                    )
                )
                for row in result.joins
            ),
        ),
    )
    atomic_write_lines(
        output_root / "mt2t_sequence_routing.tsv",
        (
            "source_ID\tdestination_ID\tstatus\treason\tevidence",
            *(
                "\t".join(
                    (
                        row.source_id,
                        row.destination_id,
                        row.status,
                        row.reason,
                        row.evidence,
                    )
                )
                for row in result.routing
            ),
        ),
    )
    atomic_write_lines(
        output_root / "mt2t_id_map.tsv",
        (
            "source_ID\tdestination_ID",
            *(f"{row.source_id}\t{row.destination_id}" for row in result.routing),
        ),
    )
    versions = dict(sorted((tool_versions or {}).items()))
    output_names = (
        "mT2T.fa",
        "retained.fa",
        "removed.fa",
        "mt2t_alignment_audit.tsv",
        "mt2t_chain_audit.tsv",
        "mt2t_overlap_edges.tsv",
        "mt2t_containment_audit.tsv",
        "mt2t_join_audit.tsv",
        "mt2t_sequence_routing.tsv",
        "mt2t_id_map.tsv",
    )
    atomic_write_lines(
        manifest_path,
        (
            "key\tvalue",
            f"phap_version\t{__version__}",
            f"fasta\t{fasta_path}",
            f"fasta_sha256\t{_sha256(fasta_path)}",
            f"paf\t{paf_path}",
            f"paf_sha256\t{_sha256(paf_path)}",
            "primary_policy\ttp:A:P_or_I",
            f"min_contig_length\t{parameters.min_contig_length}",
            f"min_alignment_block_length\t{parameters.min_alignment_block_length}",
            f"min_identity\t{parameters.min_identity:.12g}",
            f"min_shorter_coverage\t{parameters.min_shorter_coverage:.12g}",
            f"min_longer_coverage\t{parameters.min_longer_coverage:.12g}",
            f"internal_margin_ratio\t{parameters.internal_margin_ratio:.12g}",
            f"max_query_gap\t{parameters.max_query_gap}",
            f"max_target_gap\t{parameters.max_target_gap}",
            *(f"tool_{name}\t{version}" for name, version in versions.items()),
            *(
                f"output_sha256_{name}\t{_sha256(str(output_root / name))}"
                for name in output_names
            ),
        ),
    )
    return result


def _validate_routing(
    result: Mt2tAssemblyResult,
    records_by_id: Mapping[str, FastaRecord],
) -> None:
    routed = [row.source_id for row in result.routing]
    if len(routed) != len(set(routed)):
        raise AssertionError("mT2T routing contains duplicate source IDs")
    if set(routed) != set(records_by_id):
        raise AssertionError("mT2T routing does not cover every source FASTA ID")
    output_ids = {row.identifier for row in result.outputs}
    missing_destinations = sorted(
        {row.destination_id for row in result.routing} - output_ids
    )
    if missing_destinations:
        raise AssertionError(
            "mT2T routing references missing outputs: " + ", ".join(missing_destinations)
        )


def _sha256(path: str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--paf", required=True)
    parser.add_argument("--output-directory", required=True)
    parser.add_argument("--min-contig-length", type=int, default=100_000)
    parser.add_argument("--min-alignment-block-length", type=int, default=200)
    parser.add_argument("--min-identity", type=float, default=0.8)
    parser.add_argument("--min-shorter-coverage", type=float, default=0.1)
    parser.add_argument("--min-longer-coverage", type=float, default=0.05)
    parser.add_argument("--internal-margin-ratio", type=float, default=0.05)
    parser.add_argument("--max-query-gap", type=int, default=3_000_000)
    parser.add_argument("--max-target-gap", type=int, default=3_000_000)
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    parameters = Mt2tOverlapParameters(
        min_contig_length=args.min_contig_length,
        min_alignment_block_length=args.min_alignment_block_length,
        min_identity=args.min_identity,
        min_shorter_coverage=args.min_shorter_coverage,
        min_longer_coverage=args.min_longer_coverage,
        internal_margin_ratio=args.internal_margin_ratio,
        max_query_gap=args.max_query_gap,
        max_target_gap=args.max_target_gap,
    )
    run_mt2t_overlap_workflow(
        fasta_path=args.fasta,
        paf_path=args.paf,
        output_directory=args.output_directory,
        parameters=parameters,
    )


if __name__ == "__main__":
    main()
