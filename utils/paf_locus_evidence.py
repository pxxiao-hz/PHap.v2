#!/usr/bin/env python3
"""Generate strict, auditable mT2T locus evidence from complete PAF input."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple

from phap_core import __version__
from phap_core.alignment_chains import select_longest_targets
from phap_core.atomic_io import atomic_write_lines
from phap_core.fasta import read_fasta
from phap_core.locus_evidence import LocusEvidenceResult, evaluate_locus_evidence
from phap_core.locus_rescue import LocusRescueThresholds
from phap_core.paf import read_paf, select_primary_records
from phap_core.read_assignment import parse_unitig_dosages


def parse_read_support(path: str) -> Dict[str, int]:
    """Parse ``unitig_ID read_support`` with no duplicate or negative counts."""

    rows = [
        line.split()
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if not rows:
        raise ValueError("independent read-support table is empty")
    if rows[0][:2] == ["unitig_ID", "read_support"]:
        rows = rows[1:]
    support: Dict[str, int] = {}
    for row in rows:
        if len(row) != 2:
            raise ValueError("read-support rows must contain unitig_ID and read_support")
        unitig_id = row[0]
        if not unitig_id or unitig_id in support:
            raise ValueError(f"invalid or duplicate read-support unitig {unitig_id!r}")
        try:
            count = int(row[1])
        except ValueError as exc:
            raise ValueError(f"invalid read support for {unitig_id!r}") from exc
        if count < 0:
            raise ValueError(f"read support for {unitig_id!r} must be non-negative")
        support[unitig_id] = count
    return support


def parse_target_list(path: str) -> Tuple[str, ...]:
    target_rows = [
        line.strip()
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if len(target_rows) != len(set(target_rows)):
        raise ValueError("target list contains duplicate IDs")
    targets = tuple(sorted(target_rows))
    if not targets:
        raise ValueError("target list contains no IDs")
    return targets


def write_outputs(
    result: LocusEvidenceResult,
    *,
    filtered_paf: str,
    alignment_audit: str,
    candidate_audit: str,
    decision_audit: str,
    routing_audit: str,
) -> None:
    atomic_write_lines(
        filtered_paf,
        (record.to_line() for record in result.filtered_records),
    )
    atomic_write_lines(
        alignment_audit,
        [
            "source\tline\tunitig_ID\tlocus_ID\ttp\tstatus\treason\t"
            "final_status\tfinal_reason",
            *(
                "\t".join(
                    (
                        row.source,
                        str(row.line_number),
                        row.unitig_id,
                        row.locus_id,
                        row.alignment_type or ".",
                        row.status,
                        row.reason,
                        row.final_status,
                        row.final_reason,
                    )
                )
                for row in result.alignment_audit
            ),
        ],
    )
    candidate_status = {
        (decision.unitig_id, candidate.locus_id): candidate
        for decision in result.decisions
        for candidate in decision.candidates
    }
    decisions_by_unitig = {
        decision.unitig_id: decision for decision in result.decisions
    }
    candidate_lines = [
        "unitig_ID\tlocus_ID\tstrand\ttarget_start\ttarget_end\trecord_count\t"
        "matching_bases\talignment_block_bases\tidentity\tunion_query_bases\t"
        "query_coverage\tscore\tcollinear\tmeets_thresholds\t"
        "selected_for_output\tselection_reason"
    ]
    for chain in result.chains:
        candidate = candidate_status[(chain.unitig_id, chain.locus_id)]
        decision = decisions_by_unitig[chain.unitig_id]
        selected = decision.assigned_locus == chain.locus_id
        if selected:
            selection_reason = "selected_assigned_locus"
        elif decision.assigned_locus is None:
            selection_reason = "unitig_unassigned"
        else:
            selection_reason = "competing_locus_not_selected"
        candidate_lines.append(
            "\t".join(
                (
                    chain.unitig_id,
                    chain.locus_id,
                    chain.strand,
                    str(chain.target_start),
                    str(chain.target_end),
                    str(chain.record_count),
                    str(chain.matching_bases),
                    str(chain.alignment_block_bases),
                    f"{chain.identity:.12g}",
                    str(chain.union_query_bases),
                    f"{chain.query_coverage:.12g}",
                    f"{chain.score:.12g}",
                    str(chain.collinear).lower(),
                    str(candidate.meets_thresholds).lower(),
                    str(selected).lower(),
                    selection_reason,
                )
            )
        )
    atomic_write_lines(candidate_audit, candidate_lines)

    decision_lines = [
        "unitig_ID\tsource_state\tstatus\tassigned_locus\tassigned_group\t"
        "rescue_class\treason\tindependent_read_support\tbest_locus\t"
        "best_identity\tbest_query_coverage\tbest_score\tnext_best_locus\t"
        "next_best_score\tnext_best_margin"
    ]
    routing_lines = [
        "source_ID\tdestination_type\tdestination_ID\treason\tevidence"
    ]
    for decision in result.decisions:
        decision_lines.append(
            "\t".join(
                (
                    decision.unitig_id,
                    decision.source_state,
                    decision.status,
                    decision.assigned_locus or ".",
                    decision.assigned_group or ".",
                    decision.rescue_class or ".",
                    decision.reason,
                    str(decision.independent_read_support),
                    decision.best_locus or ".",
                    _optional_float(decision.best_identity),
                    _optional_float(decision.best_query_coverage),
                    _optional_float(decision.best_score),
                    decision.next_best_locus or ".",
                    _optional_float(decision.next_best_score),
                    _optional_float(decision.next_best_margin),
                )
            )
        )
        if decision.assigned_locus is not None:
            destination_type = "locus"
            destination_id = decision.assigned_locus
        else:
            destination_type = decision.status
            destination_id = decision.status
        routing_lines.append(
            "\t".join(
                (
                    decision.unitig_id,
                    destination_type,
                    destination_id,
                    decision.reason,
                    "mT2T_alignment_chain",
                )
            )
        )
    atomic_write_lines(decision_audit, decision_lines)
    atomic_write_lines(routing_audit, routing_lines)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paf", required=True, help="Complete minimap2 PAF")
    parser.add_argument("--unitig-fasta", required=True, help="Complete source unitig FASTA")
    parser.add_argument("--contig-type", required=True, help="Dosage classification table")
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument("--target-count", type=int)
    target_group.add_argument("--target-list")
    parser.add_argument("--read-support", help="Optional prevalidated unitig read support")
    parser.add_argument("--min-query-length", type=int, default=1000)
    parser.add_argument("--min-alignment-block-length", type=int, default=200)
    parser.add_argument("--max-query-gap", type=int, default=3000000)
    parser.add_argument("--max-target-gap", type=int, default=3000000)
    parser.add_argument("--min-locus-identity", type=float, default=0.8)
    parser.add_argument("--min-locus-query-coverage", type=float, default=0.05)
    parser.add_argument("--min-locus-score-margin", type=float, default=0.05)
    parser.add_argument("--min-low-coverage-read-support", type=int, default=1)
    parser.add_argument("--filtered-paf", required=True)
    parser.add_argument("--alignment-audit", required=True)
    parser.add_argument("--candidate-audit", required=True)
    parser.add_argument("--decision-audit", required=True)
    parser.add_argument("--routing-audit", required=True)
    parser.add_argument("--manifest", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    _validate_distinct_paths(args)
    _invalidate_owned_manifest(args.manifest)
    records = read_paf(args.paf)
    primary_records, _ = select_primary_records(records, require_tp=True)
    if args.target_list:
        targets = parse_target_list(args.target_list)
    else:
        if args.target_count is None:
            raise AssertionError("target count is required without a target list")
        targets = select_longest_targets(
            primary_records,
            target_count=args.target_count,
        )
    with open(args.contig_type, encoding="utf-8") as source:
        dosages, source_states = parse_unitig_dosages(source)
    fasta_unitigs = tuple(record.identifier for record in read_fasta(args.unitig_fasta))
    fasta_unitig_set = set(fasta_unitigs)
    unknown_paf_ids = sorted(
        {record.query_name for record in records} - fasta_unitig_set
    )
    if unknown_paf_ids:
        raise ValueError(
            "PAF references query IDs absent from FASTA: "
            + ", ".join(unknown_paf_ids[:5])
        )
    unknown_dosage_ids = sorted(set(dosages) - fasta_unitig_set)
    if unknown_dosage_ids:
        raise ValueError(
            "dosage table references IDs absent from FASTA: "
            + ", ".join(unknown_dosage_ids[:5])
        )
    support = parse_read_support(args.read_support) if args.read_support else {}
    unknown_support_ids = sorted(set(support) - fasta_unitig_set)
    if unknown_support_ids:
        raise ValueError(
            "read-support table references IDs absent from FASTA: "
            + ", ".join(unknown_support_ids[:5])
        )
    thresholds = LocusRescueThresholds(
        min_identity=args.min_locus_identity,
        min_query_coverage=args.min_locus_query_coverage,
        min_next_best_margin=args.min_locus_score_margin,
        min_low_coverage_read_support=args.min_low_coverage_read_support,
    )
    result = evaluate_locus_evidence(
        records,
        unitig_ids=fasta_unitigs,
        source_states=source_states,
        independent_read_support=support,
        allowed_loci=targets,
        thresholds=thresholds,
        min_query_length=args.min_query_length,
        min_alignment_block_length=args.min_alignment_block_length,
        max_query_gap=args.max_query_gap,
        max_target_gap=args.max_target_gap,
    )
    write_outputs(
        result,
        filtered_paf=args.filtered_paf,
        alignment_audit=args.alignment_audit,
        candidate_audit=args.candidate_audit,
        decision_audit=args.decision_audit,
        routing_audit=args.routing_audit,
    )
    atomic_write_lines(
        args.manifest,
        (
            "key\tvalue",
            f"phap_version\t{__version__}",
            f"paf\t{args.paf}",
            f"unitig_fasta\t{args.unitig_fasta}",
            f"contig_type\t{args.contig_type}",
            "primary_policy\ttp:A:P_or_I",
            f"allowed_loci\t{','.join(result.allowed_loci)}",
            f"min_query_length\t{args.min_query_length}",
            f"min_alignment_block_length\t{args.min_alignment_block_length}",
            f"max_query_gap\t{args.max_query_gap}",
            f"max_target_gap\t{args.max_target_gap}",
            f"min_locus_identity\t{args.min_locus_identity:.12g}",
            (
                "min_locus_query_coverage\t"
                f"{args.min_locus_query_coverage:.12g}"
            ),
            f"min_locus_score_margin\t{args.min_locus_score_margin:.12g}",
            (
                "min_low_coverage_read_support\t"
                f"{args.min_low_coverage_read_support}"
            ),
            f"read_support\t{args.read_support or '.'}",
        ),
    )


def _optional_float(value: Optional[float]) -> str:
    return "." if value is None else f"{value:.12g}"


def _validate_distinct_paths(args: argparse.Namespace) -> None:
    paths = {
        "paf": args.paf,
        "unitig_fasta": args.unitig_fasta,
        "contig_type": args.contig_type,
        "filtered_paf": args.filtered_paf,
        "alignment_audit": args.alignment_audit,
        "candidate_audit": args.candidate_audit,
        "decision_audit": args.decision_audit,
        "routing_audit": args.routing_audit,
        "manifest": args.manifest,
    }
    if args.target_list:
        paths["target_list"] = args.target_list
    if args.read_support:
        paths["read_support"] = args.read_support

    by_path: Dict[Path, str] = {}
    for role, raw_path in paths.items():
        path = Path(raw_path).resolve()
        previous_role = by_path.get(path)
        if previous_role is not None:
            raise ValueError(
                f"input and output paths must be distinct: {previous_role} and "
                f"{role} both resolve to {path}"
            )
        by_path[path] = role


def _invalidate_owned_manifest(path: str) -> None:
    manifest = Path(path)
    if not manifest.exists():
        return
    if not manifest.is_file():
        raise ValueError(f"locus manifest is not a regular file: {manifest}")
    try:
        rows = manifest.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise ValueError(f"cannot validate existing locus manifest: {manifest}") from error
    keys = {
        row.split("\t", 1)[0]
        for row in rows[1:]
        if "\t" in row
    }
    required_keys = {"phap_version", "paf", "unitig_fasta", "contig_type"}
    if not rows or rows[0] != "key\tvalue" or not required_keys.issubset(keys):
        raise ValueError(
            f"refusing to replace unrecognized locus manifest: {manifest}"
        )
    manifest.unlink()


if __name__ == "__main__":
    main()
