#!/usr/bin/env python
"""Generate an allelic-contig table from unitig-to-mT2T alignments."""

import argparse

from phap_core.allelic_bins import (
    TargetIntervalEvidence,
    build_allelic_bin_candidates,
)
from phap_core.allelic_selection import select_dosage_capacity_bins
from phap_core.atomic_io import atomic_write_lines
from phap_core.paf import parse_paf_line as parse_paf_record
from phap_core.read_assignment import parse_unitig_dosages


def parse_paf_line(
    line,
    min_align_length,
    min_unitig_length,
    source="<memory>",
    line_number=1,
):
    record = parse_paf_record(
        line,
        source=source,
        line_number=line_number,
    )

    ## filter: explicit minimap2 primary alignment
    if not record.is_primary:
        return None

    ## filter: min alignment length
    if record.alignment_block_length < min_align_length:
        return None

    ## filter: min unitig length
    if record.query_length < min_unitig_length:
        return None

    return {
        'query_id': record.query_name,
        'query_len': record.query_length,
        'query_start': record.query_start,
        'query_end': record.query_end,
        'target_id': record.target_name,
        'target_len': record.target_length,
        'target_start': record.target_start,
        'target_end': record.target_end,
        'match_len': record.matching_bases,
    }


def read_paf(file_paf, min_align_length, min_unitig_length, chr_num):
    if min_align_length < 1 or min_unitig_length < 1 or chr_num < 1:
        raise ValueError(
            "min_align_length, min_unitig_length, and chr_num must be positive"
        )
    alignments = []
    target_lengths = {}
    with open(file_paf, encoding="utf-8") as source:
        for line_number, line in enumerate(source, start=1):
            alignment = parse_paf_line(
                line,
                min_align_length,
                min_unitig_length,
                file_paf,
                line_number,
            )
            if alignment:
                alignments.append(alignment)
                previous_length = target_lengths.setdefault(
                    alignment['target_id'], alignment['target_len']
                )
                if previous_length != alignment['target_len']:
                    raise ValueError(
                        f"conflicting lengths for target {alignment['target_id']!r}"
                    )

    if len(target_lengths) > chr_num:
        raise ValueError(
            f"PAF contains {len(target_lengths)} targets, exceeding chr_num={chr_num}; "
            "provide strict locus-selected PAF evidence"
        )
    query_targets = {}
    for alignment in alignments:
        query_targets.setdefault(alignment['query_id'], set()).add(
            alignment['target_id']
        )
    ambiguous = sorted(
        query_id for query_id, targets in query_targets.items() if len(targets) > 1
    )
    if ambiguous:
        raise ValueError(
            "PAF contains unitigs assigned to multiple targets; run strict locus "
            f"selection first: {','.join(ambiguous)}"
        )
    return alignments


def build_candidates(alignments, bin_size):
    evidence = (
        TargetIntervalEvidence(
            unitig_id=alignment['query_id'],
            target_id=alignment['target_id'],
            target_length=alignment['target_len'],
            target_start=alignment['target_start'],
            target_end=alignment['target_end'],
        )
        for alignment in alignments
    )
    return build_allelic_bin_candidates(evidence, bin_size=bin_size)


def write_top_contigs_per_bin(
    file_paf,
    min_align_length,
    min_unitig_length,
    bin_size,
    ploidy,
    min_bin_support_bases,
    min_bin_coverage,
    chr_num,
    output_file,
    out_allelic_table,
    file_contig_type,
    selection_audit,
    bin_audit,
):
    with open(file_contig_type, encoding="utf-8") as dosage_source:
        dosage_by_unitig, source_states = parse_unitig_dosages(dosage_source)

    alignments = read_paf(
        file_paf,
        min_align_length,
        min_unitig_length,
        chr_num,
    )
    candidates = build_candidates(alignments, bin_size)
    decisions = select_dosage_capacity_bins(
        candidates,
        dosage_by_unitig,
        source_states,
        ploidy=ploidy,
        min_support_bases=min_bin_support_bases,
        min_bin_coverage=min_bin_coverage,
    )

    top_lines = ["#target_id\tbin_start\tbin_end\tquery_id\tmatch_len"]
    allelic_lines = []
    candidate_lines = [
        "target_id\tbin_start\tbin_end\tunitig_ID\tsource_state\tdosage\t"
        "union_support_bases\tbin_coverage\tcopy_weighted_support_bases\t"
        "status\treason\tbin_status\tbin_reason\tbest_objective\t"
        "optimal_solution_count"
    ]
    bin_lines = [
        "target_id\tbin_start\tbin_end\tploidy\tstatus\treason\t"
        "candidate_count\teligible_candidate_count\tselected_unitigs\t"
        "selected_copy_count\tremaining_capacity\tmin_support_bases\t"
        "min_bin_coverage\tobjective\tbest_objective\toptimal_solution_count"
    ]
    for decision in decisions:
        key = decision.bin_key
        selected = sorted(
            (
                row.candidate
                for row in decision.candidates
                if row.status == "selected"
            ),
            key=lambda candidate: (
                -candidate.union_support_bases,
                candidate.unitig_id,
            ),
        )
        if selected:
            for candidate in selected:
                top_lines.append(
                    f"{key.target_id}\t{key.start}\t{key.end}\t"
                    f"{candidate.unitig_id}\t{candidate.union_support_bases}"
                )
            allelic_lines.append(
                "\t".join(
                    (
                        key.target_id,
                        str(key.start),
                        str(key.end),
                        *(candidate.unitig_id for candidate in selected),
                    )
                )
            )

        for row in decision.candidates:
            candidate = row.candidate
            bin_coverage = candidate.union_support_bases / (key.end - key.start)
            candidate_lines.append(
                "\t".join(
                    (
                        key.target_id,
                        str(key.start),
                        str(key.end),
                        candidate.unitig_id,
                        row.source_state,
                        _optional_value(row.dosage),
                        str(candidate.union_support_bases),
                        f"{bin_coverage:.6f}",
                        _optional_value(row.copy_weighted_support_bases),
                        row.status,
                        row.reason,
                        decision.status,
                        decision.reason,
                        _optional_value(decision.best_objective),
                        str(decision.optimal_solution_count),
                    )
                )
            )
        eligible_count = sum(
            row.copy_weighted_support_bases is not None for row in decision.candidates
        )
        bin_lines.append(
            "\t".join(
                (
                    key.target_id,
                    str(key.start),
                    str(key.end),
                    str(ploidy),
                    decision.status,
                    decision.reason,
                    str(len(decision.candidates)),
                    str(eligible_count),
                    ",".join(decision.selected_unitigs) or ".",
                    str(decision.selected_copy_count),
                    str(decision.remaining_capacity),
                    str(decision.min_support_bases),
                    f"{decision.min_bin_coverage:.6f}",
                    "copy_weighted_target_union_bases",
                    _optional_value(decision.best_objective),
                    str(decision.optimal_solution_count),
                )
            )
        )

    atomic_write_lines(output_file, top_lines)
    atomic_write_lines(out_allelic_table, allelic_lines)
    atomic_write_lines(selection_audit, candidate_lines)
    atomic_write_lines(bin_audit, bin_lines)


def _optional_value(value):
    return "NA" if value is None else str(value)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf_file", type=str, required=True, help="PAF file")
    parser.add_argument('--min_align_length', type=int, default=1000, help="Minimum alignment length [1000]")
    parser.add_argument("--min_unitig_length", type=int, default=100000, help="Minimum unitig length [100000]")
    parser.add_argument("--bin_size", type=int, required=False, default=100000,
                        help="Bin size [100000]")
    parser.add_argument(
        "--min-bin-support-bases",
        type=int,
        default=1,
        help="Minimum target interval-union support within one bin, in bases [1]",
    )
    parser.add_argument(
        "--min-bin-coverage",
        type=float,
        default=0.0,
        help="Minimum target interval-union fraction within one bin [0.0]",
    )
    parser.add_argument(
        "--ploidy",
        "--top_n",
        dest="ploidy",
        type=int,
        default=4,
        help="Genome ploidy; --top_n is a deprecated alias [4]",
    )
    parser.add_argument('--out_top_contigs_per_bin', type=str, required=False, default='top_contigs_per_bin.txt',
                        help="Output file [top_contigs_per_bin.txt]")
    parser.add_argument('--out_allelic_table', type=str, required=False, default='allelic.ctg.table',
                        help="Output file [allelic.ctg.table]")
    parser.add_argument('--chr_num', type=int, default=12, help='number of chromosomes [12]')
    parser.add_argument('--contig_type', type=str, required=True, help='contig type based on dosage analysis')
    parser.add_argument(
        '--selection-audit',
        default='allelic_bin_candidates.tsv',
        help='Per-bin candidate routing audit [allelic_bin_candidates.tsv]',
    )
    parser.add_argument(
        '--bin-audit',
        default='allelic_bin_decisions.tsv',
        help='Per-bin capacity decision audit [allelic_bin_decisions.tsv]',
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    if (
        args.ploidy < 1
        or args.bin_size < 1
        or args.min_bin_support_bases < 1
        or not 0.0 <= args.min_bin_coverage <= 1.0
    ):
        raise SystemExit(
            'allelic_table_generate: error: invalid ploidy, bin size, or bin support threshold'
        )
    try:
        write_top_contigs_per_bin(
            args.paf_file,
            args.min_align_length,
            args.min_unitig_length,
            args.bin_size,
            args.ploidy,
            args.min_bin_support_bases,
            args.min_bin_coverage,
            args.chr_num,
            args.out_top_contigs_per_bin,
            args.out_allelic_table,
            args.contig_type,
            args.selection_audit,
            args.bin_audit,
        )
    except ValueError as error:
        raise SystemExit(f'allelic_table_generate: error: {error}') from error


if __name__ == '__main__':
    main()
