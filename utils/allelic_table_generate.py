#!/usr/bin/env python
"""Generate an allelic-contig table from unitig-to-mT2T alignments."""

import argparse
from collections import defaultdict

from phap_core.allelic_bins import (
    TargetIntervalEvidence,
    build_allelic_bin_candidates,
    select_legacy_top_n_candidates,
)
from phap_core.paf import parse_paf_line as parse_paf_record


## 解析 contig type
def parse_contig_type_based_on_dosage(file_contig_type):
    dic_dosage_contig_type = {}
    for lines in open(file_contig_type, 'r'):
        if lines.startswith('contig_ID'):
            continue
        else:
            line = lines.strip().split()
            dic_dosage_contig_type[line[0]] = line[2]
    return dic_dosage_contig_type


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
    alignments = []
    dic_scaffold_length = {}        # 存储 Scaffolds 长度
    with open(file_paf, 'r') as f:
        for line_number, line in enumerate(f, start=1):
            alignment = parse_paf_line(
                line,
                min_align_length,
                min_unitig_length,
                file_paf,
                line_number,
            )
            if alignment:                       # Skip supplementary alignment
                alignments.append(alignment)
                dic_scaffold_length[alignment['target_id']] = alignment['target_len']

    ## 获得染色体匹配数量的 scaffold id
    scffolds_chrs = sorted(dic_scaffold_length.items(), key=lambda x: x[1], reverse=True)[:chr_num]
    scffolds_chr = [scaffold_id for scaffold_id, length in scffolds_chrs]
    print(scffolds_chr)
    return alignments, scffolds_chr


def filter_longest_reference_per_contig(alignments):
    ''' 每一个 contig，只保留最长匹配的 reference '''
    contig_to_refs = defaultdict(lambda: defaultdict(int))

    # Calculate total match length per reference for each contig
    for aln in alignments:
        contig_to_refs[aln['query_id']][aln['target_id']] += aln['match_len']

    filtered_alignments = []
    for query_id, refs in contig_to_refs.items():
        # Find the reference with the longest total match length
        best_ref = max(refs, key=refs.get)
        filtered_alignments.extend(
            aln for aln in alignments if aln['query_id'] == query_id and aln['target_id'] == best_ref
        )
    return filtered_alignments


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


def group_candidates(candidates):
    grouped = defaultdict(list)
    for candidate in candidates:
        grouped[candidate.bin_key].append(candidate)
    return grouped


def write_top_contigs_per_bin(file_paf, min_align_length, min_unitig_length, bin_size, top_n, chr_num, output_file, out_allelic_table, file_contig_type):

    parse_contig_type_based_on_dosage(file_contig_type)

    alignments, scffolds_chr = read_paf(file_paf, min_align_length, min_unitig_length, chr_num)

    filtered_alignments = filter_longest_reference_per_contig(alignments)

    candidates = build_candidates(filtered_alignments, bin_size)
    selected = select_legacy_top_n_candidates(candidates, top_n=top_n)
    top_contigs_per_bin = group_candidates(selected)

    with open(output_file, 'w') as out, open(out_allelic_table, 'w') as out_allelic:
        out.write("#target_id\tbin_start\tbin_end\tquery_id\tmatch_len\n")
        for bin_key in sorted(top_contigs_per_bin):
            contigs = top_contigs_per_bin[bin_key]
            target_id = bin_key.target_id
            bin_start = bin_key.start
            bin_end = bin_key.end
            if target_id not in scffolds_chr:
                continue
            unitig_ids = []
            for candidate in contigs:
                out.write(
                    f"{target_id}\t{bin_start}\t{bin_end}\t"
                    f"{candidate.unitig_id}\t{candidate.union_support_bases}\n"
                )
                unitig_ids.append(candidate.unitig_id)
            out_allelic.write(
                '\t'.join((target_id, str(bin_start), str(bin_end), *unitig_ids))
                + '\n'
            )


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf_file", type=str, required=True, help="PAF file")
    parser.add_argument('--min_align_length', type=int, default=1000, help="Minimum alignment length [1000]")
    parser.add_argument("--min_unitig_length", type=int, default=100000, help="Minimum unitig length [100000]")
    parser.add_argument("--bin_size", type=int, required=False, default=100000,
                        help="Bin size [100000]")
    parser.add_argument("--top_n", type=int, required=False, default=4,
                        help="The maximum number of contigs corresponding to each bin."
                             " This corresponds to the number of haplotypes, for example, tetraploid is 4 [4]")
    parser.add_argument('--out_top_contigs_per_bin', type=str, required=False, default='top_contigs_per_bin.txt',
                        help="Output file [top_contigs_per_bin.txt]")
    parser.add_argument('--out_allelic_table', type=str, required=False, default='allelic.ctg.table',
                        help="Output file [allelic.ctg.table]")
    parser.add_argument('--chr_num', type=int, default=12, help='number of chromosomes [12]')
    parser.add_argument('--contig_type', type=str, required=True, help='contig type based on dosage analysis')
    args = parser.parse_args()
    return args


def main():

    args = parse_arguments()
    # 使用示例
    file_paf = args.paf_file
    min_align_length = args.min_align_length
    min_unitig_length = args.min_unitig_length
    bin_size = args.bin_size
    top_n = args.top_n
    out_top_contigs_per_bin = args.out_top_contigs_per_bin
    out_allelic_table = args.out_allelic_table
    chr_num = args.chr_num
    file_contig_type = args.contig_type

    write_top_contigs_per_bin(file_paf, min_align_length, min_unitig_length, bin_size, top_n, chr_num, out_top_contigs_per_bin, out_allelic_table, file_contig_type)


if __name__ == '__main__':
    main()
