#!/usr/bin/env python
'''
time: 2024-06-12
author: pxxiao
version: 1.0
description: 对 allelic table 重新洗牌
----\n
time: 2024-06-17
author: pxxiao
version: 2.0
description: 考虑collapsed unitig 信息，每一行 allelic table haplotype 加和不应该超过设置的 haplotype 数量。
                如果 unitig 为 diplotig，那么代表两个 haplotype，triplotig 代表三个，tetraplotig 代表四个，haplotig 代表一个。
                如果 allelic table 加起来超过 4，则需要对这一行的 allelic table 进行过滤。
                如果这个 unitig 在前面处理过，那么这个 unitig 可以过滤掉。过滤掉一个 unitig 之后，如果还是超过 4，那么再过滤一个，直到小于等于 4。
----\n
time: 2024-06-18
author: pxxiao
version: 3.0
description: 在过滤的之前，需要增加一个条件：
                如果这一行的 allelic table haplotype 加和大于4，这一行除了 new_unitig 已经在前面出现过，并且 new_unitig 在后面还会出现，
                那么这一行allelic table 的信息就是冗余的，为了避免错误，这一行可以直接过滤掉，标记为 None。
----
time: 2024-07-19
author: pxxiao
version: 3.2
description: debugs:
before:
   4731 chr08   21500000        21600000        utg000087l      utg000061l      utg000598l
   4732 chr08   21600000        21700000        utg000087l      utg000061l      utg000598l      utg002852l
   4733 chr08   21700000        21800000        utg000061l      utg000087l      utg002852l
   4734 chr08   21800000        21900000        utg000061l      utg000087l      utg001068l      utg000274l
   4735 chr08   21900000        22000000        utg000061l      utg000087l      utg001068l
   4736 chr08   22000000        22100000        utg000087l      utg001068l      utg000061l

after:
   4078 chr08   20000000        20100000        utg000061l      utg000087l      utg000274l
   4079 chr08   21700000        21800000        utg000061l      utg000274l      utg002852l
   4080 chr08   21900000        22000000        utg000061l      utg000087l      utg001068l
   4081 chr08   22000000        22100000        utg000061l      utg000087l      utg001068l
'''


import argparse
from itertools import chain

from phap_core.allelic_table import (
    AllelicTableRow,
    filter_allelic_rows,
    find_missing_bridge_unitigs,
)
from phap_core.atomic_io import atomic_write_lines
from phap_core.read_assignment import parse_unitig_dosages


### fasta to dic
def fasta_read(file_input):
    '''
    解析 FASTA 文件，存入字典 d
    :param file_input: .fasta
    :return: d: id as key, seq as value
    '''
    d = {}
    for lines in open(file_input, "r"):
        if lines.startswith(">"):
            id = lines.strip().replace(">", "")
            d[id] = []
        else:
            d[id].append(lines.strip().upper())
    for key, values in d.items():
        d[key] = "".join(values)
    return d


def parse_allelic_table(file_allelic_table):
    ''' 解析 allelic table 文件 '''
    allelic_table = []
    with open(file_allelic_table, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            unitigs = parts[3:]
            allelic_table.append((chromosome, start, end, unitigs))
    return allelic_table


## 解析 contig type
def parse_contig_type_based_on_dosage(file_contig_type):
    with open(file_contig_type) as source:
        return parse_unitig_dosages(source)


def get_unitig_lengths(unitig_ids, dic_fasta_len):
    ''' 获取 unitig 的长度信息，可以从外部文件或预定义的字典中读取 '''
    # 假设这里有一个预定义的字典存储每个 unitig 的长度
    return {unitig_id: dic_fasta_len.get(unitig_id, 0) for unitig_id in unitig_ids}


def filter_allelic_table(allelic_table, unitig_dosages, source_states, dic_fasta_len, ploidy):
    """Adapt the legacy tuple representation to the typed core implementation."""
    rows = [
        AllelicTableRow(chrom, start, end, tuple(unitigs))
        for chrom, start, end, unitigs in allelic_table
    ]
    filtered, audit = filter_allelic_rows(
        rows,
        unitig_dosages,
        source_states,
        dic_fasta_len,
        ploidy=ploidy,
    )
    legacy_rows = [
        None
        if row is None
        else (row.chromosome, row.start, row.end, list(row.unitigs))
        for row in filtered
    ]
    audit_rows = [
        (
            row.unitig_id,
            row.chromosome,
            row.start,
            row.end,
            row.source_state,
            row.reason,
        )
        for row in audit
    ]
    return legacy_rows, audit_rows


def correct_allelic_table(
    allelic_table,
    dic_fasta_len,
    dic_unitig_dosages,
    source_states,
    ploidy,
    search_range,
):
    ''' 修正 allelic table '''
    corrected_table = []
    total_bins = len(allelic_table)
    typed_rows = [
        AllelicTableRow(chromosome, start, end, tuple(unitigs))
        for chromosome, start, end, unitigs in allelic_table
    ]

    for i in range(total_bins):
        chromosome, start, end, unitigs = allelic_table[i]
        print(f"Processing bin: {chromosome}:{start}-{end}")
        original_unitigs = set(unitigs)

        # 找出当前 bin 中缺失的 unitig
        missing_unitigs = set(
            find_missing_bridge_unitigs(
                typed_rows,
                i,
                search_range=search_range,
            )
        )

        if missing_unitigs:
            print(f"Missing unitigs in bin {chromosome}:{start}-{end}: {', '.join(missing_unitigs)}")

            # 将缺失的 unitig 添加进来
            new_unitigs = original_unitigs | missing_unitigs

            print('Debugs: new_unitigs ', new_unitigs)

            corrected_table.append((chromosome, start, end, sorted(new_unitigs)))
        else:
            corrected_table.append((chromosome, start, end, sorted(original_unitigs)))

    print(corrected_table)

    return filter_allelic_table(
        corrected_table,
        dic_unitig_dosages,
        source_states,
        dic_fasta_len,
        ploidy,
    )


def write_corrected_table1(corrected_table, output_file):
    ''' 将修正后的 table 写入文件 '''
    with open(output_file, 'w') as f:
        for chromosome, start, end, unitigs in corrected_table:
            line = f"{chromosome}\t{start}\t{end}\t" + "\t".join(unitigs) + "\n"
            f.write(line)


def write_corrected_table(corrected_table, output_file):
    ''' 将修正后的 table 写入文件 '''
    lines = (
        f"{chromosome}\t{start}\t{end}\t" + "\t".join(unitigs)
        for row in corrected_table
        if row is not None
        for chromosome, start, end, unitigs in (row,)
    )
    atomic_write_lines(output_file, lines)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Refresh allelic table',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--allelic_table', type=str, help='Allelic table', required=True)
    parser.add_argument('--fasta', type=str, help='Fasta file', required=True)
    parser.add_argument('--contig_type', type=str, help='Contig type file', required=True)
    parser.add_argument('--wd', type=str, help='working directory', required=True)
    parser.add_argument(
        '--ploidy',
        '--top_n',
        dest='ploidy',
        type=int,
        help='Genome ploidy; --top_n is a deprecated alias',
        required=True,
    )
    parser.add_argument('--search_range', type=int, help='Search range', required=False, default=5)
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()
    wd = args.wd
    dic_fasta = fasta_read(args.fasta)
    dic_fasta_len = {unitig_id: len(seq) for unitig_id, seq in dic_fasta.items()}
    dic_unitig_dosages, source_states = parse_contig_type_based_on_dosage(
        args.contig_type
    )

    output_file = f'{wd}/corrected_allelic_table.txt'  # 输出文件路径

    allelic_table = parse_allelic_table(args.allelic_table)
    corrected_table, audit_rows = correct_allelic_table(
        allelic_table,
        dic_fasta_len,
        dic_unitig_dosages,
        source_states,
        args.ploidy,
        args.search_range,
    )
    write_corrected_table(corrected_table, output_file)
    audit_file = f'{wd}/allelic_table_exclusions.tsv'
    atomic_write_lines(
        audit_file,
        chain(
            ('unitig_ID\tchromosome\tstart\tend\tsource_state\treason',),
            ('\t'.join(map(str, row)) for row in sorted(audit_rows)),
        ),
    )


if __name__ == '__main__':
    main()
