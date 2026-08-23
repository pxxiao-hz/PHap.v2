#!/usr/bin/env python
'''
time: 2024-06-06
author: pxxiao
version: 1.0
description: 我们已经获得了 allelic table，dosage contig type，和 Hi-C links；利用这些信息，对 allelic unitig 分型，聚类为 4 个 group。
----
time: 2024-06-13
author: pxxiao
version: 2.0
description: 优化 Hi-C links 的分配
----
time: 2024-06-15
author: pxxiao
version: 3.0
description: 解决 allelic table 不准确的问题
            genotype 的时候，会发现 allelic table 不准确的现象，即 allelic unitig 之间存在重叠，但是却属于一个 haplotype
            针对这种情况，我们需要判断 allelic unitig 之间的信号与grouped unitig之间的信号强弱，
            如果 allelic unitig 之间的信号强于与 grouped unitig之间的信号强弱，则将allelic unitig 归为一类。
----
time: 2024-06-16
author: pxxiao
version: 4.0
description: 解决 collapsed 的问题
            目前，并没有考虑 collapsed unitig，现在需要把这个信息考虑进去。
----
time: 2024-06-18
author: pxxiao
version: 5.0
description: 解决 collapsed 的问题
            对于 allelic table 不规则的 group，即：如果有两个 group 为空，new_unitig为 triplotig，这时候不能正确分配 new_unitig。修改代码，完成这个功能。
---
version7.0:
    chr1 genotype 亲本验证结果：
        g1/log_g1_trioevale_out:N	1230	357877	0.003425
        g2/log_g2_trioevale_out:N	1418120	26237	0.018165
        g3/log_g3_trioevale_out:N	24258	219496	0.099518
        g4/log_g4_trioevale_out:N	829353	58010	0.065373
----
time: 2024-06-18
author: pxxiao
version7.0.v2
description: 进行优化
    1. 增加 flank or full 模式选择参数，如果是 flank，则 RE 酶切位点个数为两端 flank 区域的 RE 个数，而不是全长的个数。
---

To do:
    1. 有的多倍体基因组不需要 collapsed 恢复，但是其他挂载软件，例如：HapHiC 又不能解决这些基因组的挂载问题。
        这时候需要添加一个参数，跳过 collapsed 的步骤。

'''



import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import pickle
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import os

import glob

from util import *


## 解析 FASTA: ctg: [seq, seq length, RE site's count]
def parse_fasta(fasta, flank, RE='GATC'):
    """Parse input FASTA file and save sequences, lengths, and RE site counts of contigs into a dict."""
    ## 对序列计算 RE 个数
    def count_RE_sites(sequence, RE):
        ''' Count the number of restriction enzyme (RE) sites in a sequence. '''
        return sequence.count(RE)

    fa_dict = {}
    with open(fasta) as f:
        for line in f:
            if not line.strip(): continue
            if line.startswith('>'):
                ctg = line.split()[0][1:]
                fa_dict[ctg] = []
            else:
                fa_dict[ctg].append(line.strip())
    for ctg, seq_list in fa_dict.items():
        # Joining list is faster than concatenating strings
        seq = ''.join(seq_list)

        if flank is not None:
            # Ensure flank does not exceed sequence length
            flank_size = min(flank, len(seq) // 2)
            seq_segment = seq[:flank_size] + seq[-flank_size:]
        else:
            seq_segment = seq

        # Count RE sites in the determined sequence segment
        RE_sites = count_RE_sites(seq_segment, RE) + 1  # Add pseudo-count of 1 to prevent division by zero
        fa_dict[ctg] = [seq, len(seq), RE_sites]

        # Add pseudo-count of 1 to prevent division by zero (as what ALLHiC does)
        # RE_sites = count_RE_sites(seq, RE) + 1
        # fa_dict[ctg] = [seq, len(seq), RE_sites]
    return fa_dict


def adjust_hic_signals(hic_data, dic_contig_type):
    """
    Adjust Hi-C interaction signals based on the types of unitigs involved.

    Args:
    hic_data (dict): Dictionary with keys as (ctg1, ctg2) and values as interaction counts.
    unitig_types (dict): Dictionary with unitig names as keys and their types as values.

    Returns:
    dict: Adjusted Hi-C interaction data.
    """
    adjusted_hic_data = {}
    weights = {
        'haplotig': 1.0,
        'diplotig': 0.5,
        'triplotig': 0.333,
        'tetraplotig': 0.25
    }

    for (ctg1, ctg2), count in hic_data.items():
        weight1 = weights.get(dic_contig_type.get(ctg1, 'haplotig'), 1)
        weight2 = weights.get(dic_contig_type.get(ctg2, 'haplotig'), 1)
        adjusted_count = count * weight1 * weight2
        adjusted_hic_data[(ctg1, ctg2)] = adjusted_count

    return adjusted_hic_data


def parse_pickle(fa_dict, pickle_file, dic_contig_type):
    ''' 从 pickle_file 获得 full_links '''
    with open(pickle_file, 'rb') as f:
        full_link_dict = pickle.load(f)

    full_link_dict = adjust_hic_signals(full_link_dict, dic_contig_type)

    # sort by contig length
    sorted_ctg_list = sorted([(ctg, fa_dict[ctg][1]) for ctg in fa_dict], key=lambda x: x[1], reverse=True)
    # RE site
    RE_site_dict = {ctg: ctg_info[2] for ctg, ctg_info in fa_dict.items()}

    return full_link_dict, sorted_ctg_list, RE_site_dict


def load_pickle_file(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data


def get_link_list(dic_pair_hic):
    ''' 将 contig pair 的 Hi-C 信号，转换为 contig 的 Hi-C 信号 '''
    dic_contig_hic = defaultdict(list)
    ## 处理 contig pair links dic
    for (contig1, contig2), links in dic_pair_hic.items():
        dic_contig_hic[contig1].append((contig2, links))
        if contig1 != contig2:
            dic_contig_hic[contig2].append((contig1, links))

    # 对每个 contig 的 links 按照从大到小排序
    for contig, links in dic_contig_hic.items():
        links.sort(key=lambda x: x[1], reverse=True)

    return dic_contig_hic


def parse_contig_type(file_contig_type):
    df = pd.read_csv(file_contig_type, delim_whitespace=True)
    dic_contig_type = dict(zip(df['contig_ID'], df['contig_type']))
    return dic_contig_type


def unitig_overlap_ratio(file_allelic_table):
    ''' 根据 allelic table 获得 unitig 之间的 overlap ratio '''
    def read_allelic_table(file_path):
        '''
        读取 allelic table 文件并解析数据，返回每个 unitig 的 bin 覆盖情况。
        :param file_path:
        :return: 包含每个 unitig 覆盖 bin 的字典
        '''
        unitig_bins = {}

        with open(file_path, 'r') as file:
            for line in file:
                columns = line.strip().split()
                chrom = columns[0]
                start = int(columns[1])
                end = int(columns[2])
                unitigs = columns[3:]

                for unitig in unitigs:
                    if unitig not in unitig_bins:
                        unitig_bins[unitig] = set()
                    unitig_bins[unitig].add((start, end))

        return unitig_bins

    def calculate_overlap_ratio(unitig_bins, unitig1, unitig2):
        '''
        计算两个 unitig 之间的 overlap ratio
        :param unitig_bins: 包含每个 unitig 覆盖 bin 的字典
        :param unitig1: 第一个 unitig
        :param unitig2: 第二个 unitig
        :return: 两个 overlap ratio，(ratio1, ratio2)
        '''
        bins1 = unitig_bins.get(unitig1, set())
        bins2 = unitig_bins.get(unitig2, set())

        if not bins1 or not bins2:
            return 0, 0

        overlap_bins = bins1 & bins2
        overlap_count = len(overlap_bins)

        ratio1 = overlap_count / len(bins1) if bins1 else 0
        ratio2 = overlap_count / len(bins2) if bins2 else 0

        return ratio1, ratio2

    unitig_bins = read_allelic_table(file_allelic_table)

    overlap_ratios = {}
    unitigs = list(unitig_bins.keys())

    for i in range(len(unitigs)):
        for j in range(i + 1, len(unitigs)):
            unitig1 = unitigs[i]
            unitig2 = unitigs[j]
            ratio1, ratio2 = calculate_overlap_ratio(unitig_bins, unitig1, unitig2)
            overlap_ratios[(unitig1, unitig2)] = ratio1
            overlap_ratios[(unitig2, unitig1)] = ratio2

    return overlap_ratios


def genotype_allelic_table(file_allelic_table, dic_contig_type, dic_contig_hic, RE_site_dict):
    def get_contig_types(contigs, dic_contig_type):
        return [dic_contig_type.get(contig, 'unknown') for contig in contigs]

    def find_group_for_new_unitig1(new_unitig, groups, dic_contig_hic, dic_contig_type, dic_RE_site):
        ''' 根据 Hi-C links 的密度找到适合的 group 分配 new_unitig '''

        ## 获取 new_unitig 的类型
        contig_type = dic_contig_type.get(new_unitig, 'haplotig')
        # 所需要的 group 数量
        required_group_count = {'haplotig': 1, 'diplotig': 2, 'triplotig': 3, 'tetraplotig': 4}.get(contig_type, 1)

        ## 计算每个 group 的 Hi-C links 总数和密度
        group_density = []
        for idx, group in enumerate(groups):
            total_hic_links = 0
            total_re_sites = 0

            # 计算该 group 的 Hi-C links 总数
            for unitig in group:
                hic_links = dic_contig_hic.get(unitig, [])
                hic_count = sum(count for u, count in hic_links if u == new_unitig)
                total_hic_links += hic_count
                total_re_sites += dic_RE_site.get(unitig, 1)  # 默认每个 unitig 至少有一个酶切位点，避免除零错误

            # 计算 Hi-C links 密度
            if total_hic_links > 0:
                density = total_hic_links / total_re_sites
            else:
                density = 0

            group_density.append((idx, density))

        # 筛选出密度最高的 groups
        group_density.sort(key=lambda x: x[1], reverse=True)

        selected_groups = [
            idx for idx, density in group_density if density > 0
        ][:required_group_count]

        # Collapsed unitigs require enough independently supported groups.
        if len(selected_groups) == required_group_count:
            return selected_groups
        else:
            return None

    def find_group_for_new_unitig(new_unitig, groups, dic_contig_hic, dic_contig_type, threshold=50):

        ## new_unitig 的类型
        contig_type = dic_contig_type.get(new_unitig, 'haplotig')

        ## 找到每个 group 中最后一个 unitig
        # last_unitigs = [group[-1] if group else None for group in groups]

        valid_last_unitigs = []
        for group in groups:
            last_unitig = group[-1] if group else None
            # 检查并处理 collapsed unitig
            while last_unitig and dic_contig_hic.get(last_unitig) in  ['diplotig', 'triplotig', 'tetraplotig']:
                # 如果当前 last_unitig 是 collapsed unitig，则回溯找前一个 unitig
                last_index = group.index(last_unitig)
                last_unitig = group[last_index - 1] if last_index > 0 else None
            valid_last_unitigs.append(last_unitig)

        ## 计算 new_unitig 与 grouped_last_unitig 的 Hi-C links
        last_link_counts = []
        for last_unitig in valid_last_unitigs:
            if last_unitig:
                links = dic_contig_hic.get(last_unitig, [])     # last unitig 的 Hi-C links 列表
                link_count = next((count for u, count in links if u == new_unitig), 0)      # last unitig 与 new_unitig 的 Hi-C link
                last_link_counts.append(link_count)
            else:
                last_link_counts.append(0)

        ## 计算 new_unitig 与所有 group 的 link 数量
        total_links_per_group = {}
        # 遍历所有的 groups
        for idx, group in enumerate(groups):
            total_links = 0
            # 遍历 group 中的每个 unitig
            for unitig in group:
                # 获取 unitig 的 Hi-C 链接数据
                links = dic_contig_hic.get(unitig, [])
                # 累加与 new_unitig 的 links
                total_links += sum(count for u, count in links if u == new_unitig)
            # 将计算结果存储到字典中
            total_links_per_group[idx] = total_links

        ## 计算 new_unitig 与所有 group 的 links density


        # debugs:
        if new_unitig == 'utg000961l':
            print(f'Debugs contig_type {contig_type}')
            print(f'Debugs xpx last_link_counts {last_link_counts}')    # [14, 0, 0]

        ## 如果 group 为空，并且满足 collapsed unitig 对应的数量，则直接返回 empty group 的索引；否则，返回 empty group 的索引和 Hi-C 匹配的 group 索引
        empty_groups = [i for i, group in enumerate(groups) if not group]
        print(f'Debugs xpx empty_groups {empty_groups}')
        required_empty_count = {
            'haplotig': 1,
            'diplotig': 2,
            'triplotig': 3,
            'tetraplotig': 4
        }.get(contig_type, 1)

        result_indices = []

        if len(empty_groups) >= required_empty_count:
            return empty_groups[:required_empty_count]      # tetraplotig 的时候，是不是直接返回四个group？
        else:
            result_indices.extend(empty_groups)

        print(f'Debugs result_indices-1 {new_unitig} {result_indices}')

        ## 根据 unitig 类型确定需要的 group 数量
        required_group_count = {
            'haplotig': 1,
            'diplotig': 2,
            'triplotig': 3,
            'tetraplotig': 4
        }.get(contig_type, 1)

        ## 检查 new_unitig 与 group 的所有 Hi-C 数量差异
        ## 找到链接数最多的group
        sorted_groups = sorted(range(len(total_links_per_group)), key=lambda i: total_links_per_group[i], reverse=True)
        best_group_indices = sorted_groups[:required_group_count]

        ## 如果 group 为空，并且满足 collapsed unitig 对应的数量，则直接返回 empty group 的索引；否则，返回 empty group 的索引和 Hi-C 匹配的 group 索引
        empty_groups = [i for i, group in enumerate(groups) if not group]
        if len(empty_groups) >= required_group_count:
            return empty_groups[:required_group_count]

        ## 如果最高链接数大于阈值，返回链接数最高的group，否则返回None
        if all(total_links_per_group[idx] > threshold for idx in best_group_indices):
            return best_group_indices
        else:
            print(
                f"Debug: No Hi-C signal found for {new_unitig} meeting the threshold with any group. Not adding to any group.")
            return None

    def get_remaining_to_original_index(groups, processed_groups):
        """
        为 remaining_groups 中的每个组生成唯一标识符，并映射到原始组的索引。
        """
        remaining_groups = [groups[i] for i in range(len(groups)) if i not in processed_groups]
        remaining_to_original_index = {}
        used_indices = set()

        for i, group in enumerate(remaining_groups):
            for j, g in enumerate(groups):
                if g == group and j not in used_indices:
                    remaining_to_original_index[i] = j
                    used_indices.add(j)
                    break  # 只匹配第一个未使用的组

        return remaining_groups, remaining_to_original_index

    g1, g2, g3, g4 = [], [], [], []
    groups = [g1, g2, g3, g4]
    dosage = {
        'haplotig': 1,
        'diplotig': 2,
        'triplotig': 3,
        'tetraplotig': 4,
    }
    first_line = True
    with open(file_allelic_table) as allelic_table:
        for lines in allelic_table:
            line = lines.strip().split()
            if not line:
                continue
            unitigs = line[3:]
            unitig_types = get_contig_types(unitigs, dic_contig_type)

            if first_line:
                total_dosage = sum(dosage.get(unitig_type, 0) for unitig_type in unitig_types)
                if total_dosage > 4 or any(unitig_type not in dosage for unitig_type in unitig_types):
                    raise ValueError(
                        f'Invalid dosage composition in allelic table row: {lines.rstrip()}'
                    )
                group_index = 0
                for unitig, unitig_type in zip(unitigs, unitig_types):
                    for _ in range(dosage[unitig_type]):
                        groups[group_index].append(unitig)
                        group_index += 1
                first_line = False
                continue

            processed_groups = set()
            processed_unitigs = set()
            for unitig in unitigs:
                for group_index, group in enumerate(groups):
                    if unitig in group:
                        processed_unitigs.add(unitig)
                        processed_groups.add(group_index)

            for unitig in unitigs:
                if unitig in processed_unitigs:
                    continue
                remaining_groups, remaining_to_original_index = get_remaining_to_original_index(
                    groups, processed_groups
                )
                target_group_indices = find_group_for_new_unitig1(
                    unitig,
                    remaining_groups,
                    dic_contig_hic,
                    dic_contig_type,
                    RE_site_dict,
                )

                if target_group_indices is not None:
                    for target_group_index in target_group_indices:
                        original_group_index = remaining_to_original_index[target_group_index]
                        if unitig not in groups[original_group_index]:
                            groups[original_group_index].append(unitig)
                        processed_groups.add(original_group_index)
                    processed_unitigs.add(unitig)
                elif not all(groups):
                    # Projection evidence may initialize empty chromosome-end
                    # groups even when the unitig has no observed Hi-C contact.
                    required = dosage.get(dic_contig_type.get(unitig), 1)
                    empty_groups = [group for group in groups if not group]
                    for group in empty_groups[:required]:
                        group.append(unitig)

    return g1, g2, g3, g4


def list_to_file(list_1, file_name):
    out = open(file_name, 'w')
    list_1 = sorted(set(list_1))
    for i in list_1:
        out.write(i + '\n')
    out.close()


def list_to_cluster_file(list_1, list_name, file='group.cluster.txt'):
    out = open(file, 'a+')
    list_1 = sorted(set(list_1))
    out.write(list_name + '\t' + str(len(list_1)) + '\t' + ' '.join(list_1) + '\n')
    out.close()


def write_combined_genotypes_to_file(g1, g2, g3, g4, filename):
    """ 将 g1, g2, g3, g4 写入同一个文件，每一列是一个 genotype """
    # 找到最长的 list 的长度
    max_length = max(len(g1), len(g2), len(g3), len(g4))

    # 填充较短的 list 使其长度相同
    g1.extend([''] * (max_length - len(g1)))
    g2.extend([''] * (max_length - len(g2)))
    g3.extend([''] * (max_length - len(g3)))
    g4.extend([''] * (max_length - len(g4)))

    # 将组合后的内容写入文件
    with open(filename, 'w') as f:
        f.write('g1' + '\t' + 'g2' + '\t' + 'g3' + '\t' + 'g4' + '\n')
        for row in zip(g1, g2, g3, g4):
            f.write('\t'.join(row) + '\n')


def output_fa_from_list(list_1, dic_fasta, file_name):
    out = open(file_name, 'w')
    list_1 = sorted(set(list_1))
    for i in list_1:
        if i in dic_fasta:
            out.write('>' + i + '\n' + dic_fasta[i] + '\n')
        else:
            print(f"Warning: {i} not found in the fasta dictionary.")
    # for i in list_1:
    #     out.write('>' + i + '\n' + dic_fasta[i] + '\n')
    out.close()


def rename_id(g_fasta, merge):
    '''
    Mark 单倍型 ID，添加到 merge.fa 后面
    例如：g1.fa，ID 更新为>g1_unitigID
    '''
    bn = os.path.basename(g_fasta).replace('.fa', '')
    out = open(merge, 'a+')
    for lines in open(g_fasta, 'r'):
        if lines.startswith('>'):
            out.write('>' + bn + '_' + lines.strip().replace('>', '') + '\n')
        else:
            out.write(lines)


def false_genotype_unitig(g_list, dic_fasta, dic_contig_hic):
    ''' 被错误分型的 unitig '''
    corrected_list = []
    for u in g_list:
        if u not in dic_contig_hic:
            print(f'Debugs: No hi-c link, false genotype {u}')
        else:
            corrected_list.append(u)
    return corrected_list


def cluster(
    fasta,
    full_links,
    flank,
    contig_type,
    allelic_table_file,
    wd,
    ploidy=4,
    balance_weight=1.0,
    max_refinement_rounds=10,
    max_phase_block_rounds=10,
    min_phase_block_gain=0.005,
    max_phase_boundary_relaxations=1,
    max_phase_boundary_overlap=0.30,
    max_constraint_relaxation_rounds=50,
    min_constraint_relaxation_gain=0.0,
    min_constraint_relaxation_links=5.0,
    min_constraint_relaxation_margin=0.10,
    max_constraint_relaxations_per_move=2,
    max_constraint_relaxation_overlap=0.30,
    max_phase_interval_rounds=4,
    max_phase_interval_relaxations=2,
    max_backtracks=1_000_000,
    constraint_relaxation='weakest',
    protected_long_unitig_length=5_000_000,
    protected_length_ratio=5.0,
    protected_short_overlap=0.50,
    hic_link_normalization='dosage',
):
    cluster_script = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        '..',
        'utils',
        'cluster_allelic_unitigs_v2.py',
    )
    command = [
        sys.executable,
        cluster_script,
        '--fasta', fasta,
        '--full-links', full_links,
        '--contig-type', contig_type,
        '--allelic-table', allelic_table_file,
        '--output-dir', wd,
        '--ploidy', str(ploidy),
        '--balance-weight', str(balance_weight),
        '--max-refinement-rounds', str(max_refinement_rounds),
        '--max-phase-block-rounds', str(max_phase_block_rounds),
        '--min-phase-block-gain', str(min_phase_block_gain),
        '--max-phase-boundary-relaxations', str(max_phase_boundary_relaxations),
        '--max-phase-boundary-overlap', str(max_phase_boundary_overlap),
        '--max-constraint-relaxation-rounds', str(max_constraint_relaxation_rounds),
        '--min-constraint-relaxation-gain', str(min_constraint_relaxation_gain),
        '--min-constraint-relaxation-links', str(min_constraint_relaxation_links),
        '--min-constraint-relaxation-margin', str(min_constraint_relaxation_margin),
        '--max-constraint-relaxations-per-move', str(max_constraint_relaxations_per_move),
        '--max-constraint-relaxation-overlap', str(max_constraint_relaxation_overlap),
        '--max-phase-interval-rounds', str(max_phase_interval_rounds),
        '--max-phase-interval-relaxations', str(max_phase_interval_relaxations),
        '--max-backtracks', str(max_backtracks),
        '--constraint-relaxation', constraint_relaxation,
        '--protected-long-unitig-length', str(protected_long_unitig_length),
        '--protected-length-ratio', str(protected_length_ratio),
        '--protected-short-overlap', str(protected_short_overlap),
        '--hic-link-normalization', hic_link_normalization,
    ]
    if flank is not None:
        command.extend(['--flank', str(flank)])
    subprocess.run(command, check=True)


def run_logged_command(command, stdout_path, stderr_path):
    with open(stdout_path, 'w') as stdout, open(stderr_path, 'w') as stderr:
        subprocess.run(command, stdout=stdout, stderr=stderr, check=True)


def parse_arguments():
    parser = argparse.ArgumentParser("Haplotype clustering of autopolyploid genome.")

    parser.add_argument('--p_utg', required=True, type=str, help='Path to p_utg file')
    parser.add_argument('--mT2T', required=True, type=str, help='Path to mT2T file or reference genome')
    parser.add_argument('--contig_type', required=True, type=str, help='Path to contig type from dosage analysis')
    parser.add_argument('--threads', type=int, default=10, help='The number of threads [10]')
    parser.add_argument('--paf', type=str, help='Reuse an existing raw p_utg vs mT2T PAF')
    parser.add_argument('--gfa', type=str, help='hifiasm GFA for graph-aware conflict filtering')
    parser.add_argument('--output_dir', type=str, default='02.cluster.v2', help='Output directory [02.cluster.v2]')
    parser.add_argument('--table_only', action='store_true', help='Stop after generating the v2 allelic table and QC files')
    parser.add_argument(
        '--stop_after', choices=['table', 'chr_seq', 'cluster', 'recluster', 'rescue'],
        default='rescue', help='Stop after the selected completed workflow stage [rescue]'
    )

    find_longest = parser.add_argument_group('>>> Parameters for find longest subsequences')
    find_longest.add_argument('--min_align_length', type=int, default=1000, help="Minimum alignment length [1000]")
    find_longest.add_argument("--min_unitig_length", type=int, default=20000, help="Deprecated v1 option; use --min_query_length")
    find_longest.add_argument('--min_alignment_distance', type=int, default=500000, help='The minimum distance between two alignments [500000]')
    find_longest.add_argument('--min_match_ratio', type=float, default=0.05, help='Minimum match ratio [0.05]')
    find_longest.add_argument('--min_lis_size', type=int, default=5, help='Minimum number of alignments in a LIS to be considered for mergeing [5]')
    find_longest.add_argument('--min_lis_length', type=int, default=1000000, help='Minimum length of alignments in a LIS to be considered for mergeing [1000000]')
    find_longest.add_argument('--max_lis_distance', type=int, default=3000000, help='Maximum distance between LIS to be considered for merge [3000000]')

    allelic_table = parser.add_argument_group('>>> Allelic table generated')
    allelic_table.add_argument("--bin_size", type=int, default=100000, help="Bin size [100000]")
    allelic_table.add_argument('--chr_num', type=int, default=12, help='The number of chromosomes [12]')
    allelic_table.add_argument('--top_n', type=int, default=4, help='The number of haplotypes, for example, tetraploid is 4 [4]')
    allelic_table.add_argument('--search_range', type=int, default=20, help='Search range [20]')
    allelic_table.add_argument('--min_identity', type=float, default=0.90)
    allelic_table.add_argument('--min_alignment_mapq', type=int, default=20)
    allelic_table.add_argument('--min_query_length', type=int, default=20000)
    allelic_table.add_argument('--min_chain_aligned_bp', type=int, default=20000)
    allelic_table.add_argument('--min_query_coverage', type=float, default=0.30)
    allelic_table.add_argument('--min_target_span_coverage', type=float, default=0.30)
    allelic_table.add_argument('--min_sparse_query_span_coverage', type=float, default=0.80)
    allelic_table.add_argument('--min_sparse_aligned_bp', type=int, default=1000000)
    allelic_table.add_argument('--min_sparse_anchors', type=int, default=20)
    allelic_table.add_argument('--min_sparse_mean_mapq', type=float, default=20.0)
    allelic_table.add_argument('--min_fragmented_aligned_bp', type=int, default=500000)
    allelic_table.add_argument('--min_fragmented_anchors', type=int, default=20)
    allelic_table.add_argument('--min_fragmented_mean_mapq', type=float, default=20.0)
    allelic_table.add_argument('--min_fragmented_reference_margin', type=float, default=0.20)
    allelic_table.add_argument('--min_chain_identity', type=float, default=0.90)
    allelic_table.add_argument('--min_reference_margin', type=float, default=0.05)
    allelic_table.add_argument('--max_overlap', type=int, default=10000)
    allelic_table.add_argument('--max_projection_gap', type=int, default=20000)
    allelic_table.add_argument(
        '--max_projection_blocks', type=int, default=10,
        help='QC threshold for fragmented unitigs; does not reject the unitig [10]'
    )
    allelic_table.add_argument('--min_projection_aligned_bp', type=int, default=10000)
    allelic_table.add_argument('--min_projection_coverage', type=float, default=0.70)
    allelic_table.add_argument('--min_block_identity', type=float, default=0.90)
    allelic_table.add_argument('--min_block_mapq', type=float, default=20.0)
    allelic_table.add_argument('--min_block_collinearity', type=float, default=0.80)
    allelic_table.add_argument('--min_anchor_block_aligned_bp', type=int, default=100000)
    allelic_table.add_argument('--min_anchor_query_fraction', type=float, default=0.005)
    allelic_table.add_argument('--min_segment_length', type=int, default=10000)
    allelic_table.add_argument('--min_path_envelope_query_length', type=int, default=5000000)
    allelic_table.add_argument('--min_path_envelope_query_coverage', type=float, default=0.25)
    allelic_table.add_argument('--min_path_envelope_query_span_coverage', type=float, default=0.80)
    allelic_table.add_argument('--min_path_envelope_target_coverage', type=float, default=0.20)
    allelic_table.add_argument('--max_path_envelope_span_ratio', type=float, default=2.0)
    allelic_table.add_argument('--min_long_path_pair_overlap', type=int, default=10000)
    allelic_table.add_argument('--min_over_capacity_pair_overlap', type=int, default=10000)
    allelic_table.add_argument('--min_over_capacity_pair_short_coverage', type=float, default=0.50)
    allelic_table.add_argument('--min_deferred_over_capacity_bp', type=int, default=100000)
    allelic_table.add_argument('--max_deferred_over_capacity_query_length', type=int, default=5000000)
    allelic_table.add_argument('--max_deferred_over_capacity_preferred_fraction', type=float, default=0.0)
    allelic_table.add_argument('--length_prior_scale', type=float, default=1000000.0)
    allelic_table.add_argument('--unknown_dosage_policy', choices=['exclude', 'haplotig'], default='exclude')
    allelic_table.add_argument('--over_capacity_policy', choices=['omit', 'confident', 'best'], default='confident')
    allelic_table.add_argument('--min_resolution_margin', type=float, default=0.10)

    chromosome_assignment = parser.add_argument_group('>>> Chromosome assignment')
    chromosome_assignment.add_argument(
        '--min_chromosome_margin', type=float, default=0.05,
        help='Minimum best-vs-second chain score margin for chromosome rescue [0.05]'
    )

    cluster = parser.add_argument_group('>>> Cluster based on Hi-C links')
    cluster.add_argument('--full_links', type=str, help='Path to full links pickle')
    cluster.add_argument('--flank', type=int, help='Flank')
    cluster.add_argument(
        '--hic_link_normalization', choices=['dosage', 'raw'], default='dosage',
        help=(
            'Use dosage-normalized Hi-C counts in stages 03-05, or raw read-pair '
            'counts [dosage]'
        )
    )
    cluster.add_argument(
        '--cluster_balance_weight', type=float, default=1.0,
        help='Weight of normalized group-bp imbalance relative to Hi-C cohesion [1.0]'
    )
    cluster.add_argument(
        '--cluster_refinement_rounds', type=int, default=10,
        help='Maximum rounds for single-unit and Kempe-component refinement [10]'
    )
    cluster.add_argument('--cluster_phase_block_rounds', type=int, default=10)
    cluster.add_argument('--cluster_min_phase_block_gain', type=float, default=0.005)
    cluster.add_argument('--cluster_max_phase_boundary_relaxations', type=int, default=1)
    cluster.add_argument('--cluster_max_phase_boundary_overlap', type=float, default=0.30)
    cluster.add_argument('--cluster_max_constraint_relaxation_rounds', type=int, default=50)
    cluster.add_argument('--cluster_min_constraint_relaxation_gain', type=float, default=0.0)
    cluster.add_argument('--cluster_min_constraint_relaxation_links', type=float, default=5.0)
    cluster.add_argument('--cluster_min_constraint_relaxation_margin', type=float, default=0.10)
    cluster.add_argument('--cluster_max_constraint_relaxations_per_move', type=int, default=2)
    cluster.add_argument('--cluster_max_constraint_relaxation_overlap', type=float, default=0.30)
    cluster.add_argument('--cluster_max_phase_interval_rounds', type=int, default=4)
    cluster.add_argument('--cluster_max_phase_interval_relaxations', type=int, default=2)
    cluster.add_argument(
        '--cluster_max_backtracks', type=int, default=1000000,
        help='Maximum exact constraint-search backtracks [1000000]'
    )
    cluster.add_argument(
        '--cluster_constraint_relaxation',
        choices=['weakest', 'fail'],
        default='weakest',
        help='Resolve globally impossible allelic constraints by relaxing the weakest edge, or fail [weakest]'
    )
    cluster.add_argument('--cluster_protected_long_unitig_length', type=int, default=5000000)
    cluster.add_argument('--cluster_protected_length_ratio', type=float, default=5.0)
    cluster.add_argument('--cluster_protected_short_overlap', type=float, default=0.50)

    recluster = parser.add_argument_group('>>> Recluster based on Hi-C links')
    recluster.add_argument(
        '--recluster_min_adjusted_links', type=float, default=5.0,
        help='Minimum Hi-C support on the selected normalization scale [5.0]'
    )
    recluster.add_argument(
        '--recluster_min_group_margin', type=float, default=0.10,
        help='Minimum weakest-selected vs strongest-unselected density margin [0.10]'
    )
    recluster.add_argument('--recluster_min_allelic_block_anchors', type=int, default=1)
    recluster.add_argument(
        '--recluster_max_allelic_block_configurations', type=int, default=256
    )
    recluster.add_argument(
        '--recluster_max_allelic_block_movable_length', type=int, default=5000000,
        help=(
            'Protect unitigs at least this long from local allelic-block moves; '
            'zero disables the protection [5000000]'
        )
    )
    recluster.add_argument(
        '--recluster_min_group_bp_ratio', type=float, default=0.25,
        help=(
            'Reject catastrophic haplotype-group collapse below this fraction '
            'of median group bp; zero disables [0.25]'
        )
    )
    recluster.add_argument(
        '--recluster_min_assigned_fraction', type=float, default=0.0,
        help='Optional minimum fraction of group-linked Hi-C assigned to selected groups [0.0]'
    )
    recluster.add_argument(
        '--recluster_max_rounds', type=int, default=10,
        help='Maximum high-confidence propagation rounds [10]'
    )
    recluster.add_argument(
        '--recluster_refinement_rounds', type=int, default=4,
        help='Synchronous stability checks corresponding to recluster runs 2-5 [4]'
    )
    recluster.add_argument(
        '--recluster_low_confidence_policy', choices=['defer', 'best'], default='defer',
        help='Keep ambiguous unitigs unassigned or use their best available groups [defer]'
    )
    recluster.add_argument(
        '--recluster_unknown_dosage_policy', choices=['error', 'haplotig'], default='haplotig',
        help='Handle chromosome unitigs absent from the dosage table [haplotig]'
    )
    recluster.add_argument(
        '--recluster_seed_review', choices=['weak', 'off'], default='weak',
        help='Re-evaluate weak 03.cluster seeds with chromosome Hi-C, or keep every seed fixed [weak]'
    )
    recluster.add_argument(
        '--recluster_trusted_seed_bases', default='hic_supported',
        help='Comma-separated 03.cluster assignment bases retained as immutable anchors [hic_supported]'
    )
    recluster.add_argument(
        '--recluster_reviewed_seed_fallback', choices=['retain', 'defer'], default='retain',
        help='Retain an original seed group when Hi-C review is inconclusive, or defer it [retain]'
    )

    rescue = parser.add_argument_group('>>> Rescue chromosome-unassigned unitigs')
    rescue.add_argument(
        '--rescue_min_adjusted_links', type=float, default=5.0,
        help='Minimum support on the selected Hi-C normalization scale [5.0]'
    )
    rescue.add_argument(
        '--rescue_min_chromosome_margin', type=float, default=0.10,
        help='Minimum normalized best-vs-second chromosome margin [0.10]'
    )
    rescue.add_argument(
        '--rescue_min_group_margin', type=float, default=0.10,
        help='Minimum weakest-selected vs strongest-unselected group margin [0.10]'
    )
    rescue.add_argument(
        '--rescue_min_assigned_fraction', type=float, default=0.0,
        help='Optional minimum fraction of group-linked Hi-C in selected groups [0.0]'
    )
    rescue.add_argument(
        '--rescue_max_rounds', type=int, default=10,
        help='Maximum high-confidence rescue propagation rounds [10]'
    )
    rescue.add_argument(
        '--rescue_low_confidence_policy', choices=['defer', 'best'], default='defer',
        help='Keep ambiguous candidates unassigned or use best available groups [defer]'
    )
    rescue.add_argument(
        '--rescue_unknown_dosage_policy',
        choices=['defer', 'error', 'haplotig'], default='defer',
        help='Defer unsupported dosage types, stop, or treat them as haplotigs [defer]'
    )

    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    script_realpath = os.path.dirname(os.path.realpath(__file__))
    utils_realpath = os.path.join(script_realpath, '..', 'utils')

    cwd = os.getcwd()       # current working dir
    for path_argument in (
        'p_utg', 'mT2T', 'contig_type', 'paf', 'gfa', 'full_links'
    ):
        value = getattr(args, path_argument, None)
        if value:
            setattr(args, path_argument, os.path.abspath(value))
    output_root = os.path.abspath(args.output_dir)
    stop_after = 'table' if args.table_only else args.stop_after
    stage_order = {
        'table': 1, 'chr_seq': 2, 'cluster': 3, 'recluster': 4, 'rescue': 5
    }
    if stage_order[stop_after] >= stage_order['cluster'] and not args.full_links:
        raise SystemExit('--full_links is required when running through cluster')
    ### Step 1: p_utg vs mT2T
    step1_dir = os.path.join(output_root, '01.putg_vs_mT2T')
    os.makedirs(step1_dir, exist_ok=True)

    paf_file = os.path.abspath(args.paf) if args.paf else os.path.join(step1_dir, 'p_utg_vs_mT2T.paf')
    best_paf_file = os.path.join(step1_dir, 'putg_vs_mT2T.collinear.paf')
    chain_qc_file = os.path.join(step1_dir, 'collinear_chain.qc.tsv')
    selection_summary_file = os.path.join(step1_dir, 'chromosome_selection.summary.json')
    allelic_table_file = os.path.join(step1_dir, 'corrected_allelic_table.txt')

    try:
        # Minimap2 alignment
        if not args.paf and not os.path.exists(paf_file):
            with open(paf_file, 'w') as paf_output:
                subprocess.run([
                    'minimap2', '-cx', 'asm5', '-t', str(args.threads),
                    args.mT2T, args.p_utg,
                ], stdout=paf_output, check=True)

        subprocess.run([
            sys.executable, os.path.join(utils_realpath, 'find_collinear_chains.py'),
            '--paf', paf_file, '--output', best_paf_file, '--qc', chain_qc_file,
            '--summary', selection_summary_file,
            '--min-alignment-length', str(args.min_align_length),
            '--min-identity', str(args.min_identity),
            '--min-alignment-mapq', str(args.min_alignment_mapq),
            '--min-query-length', str(args.min_query_length),
            '--min-chain-aligned-bp', str(args.min_chain_aligned_bp),
            '--min-query-coverage', str(args.min_query_coverage),
            '--min-target-span-coverage', str(args.min_target_span_coverage),
            '--min-sparse-query-span-coverage', str(args.min_sparse_query_span_coverage),
            '--min-sparse-aligned-bp', str(args.min_sparse_aligned_bp),
            '--min-sparse-anchors', str(args.min_sparse_anchors),
            '--min-sparse-mean-mapq', str(args.min_sparse_mean_mapq),
            '--min-fragmented-aligned-bp', str(args.min_fragmented_aligned_bp),
            '--min-fragmented-anchors', str(args.min_fragmented_anchors),
            '--min-fragmented-mean-mapq', str(args.min_fragmented_mean_mapq),
            '--min-fragmented-reference-margin', str(args.min_fragmented_reference_margin),
            '--min-chain-identity', str(args.min_chain_identity),
            '--min-reference-margin', str(args.min_reference_margin),
            '--max-overlap', str(args.max_overlap),
            '--max-gap', str(args.max_lis_distance)
        ], check=True)

        table_command = [
            sys.executable, os.path.join(utils_realpath, 'allelic_table_generate_v2.py'),
            '--paf', best_paf_file, '--contig-type', args.contig_type,
            '--output', allelic_table_file,
            '--projections', os.path.join(step1_dir, 'unitig_projections.tsv'),
            '--qc', os.path.join(step1_dir, 'allelic_table.qc.tsv'),
            '--rejected', os.path.join(step1_dir, 'rejected_projections.tsv'),
            '--pairs', os.path.join(step1_dir, 'allelic_pairs.tsv'),
            '--summary', os.path.join(step1_dir, 'allelic_table.summary.json'),
            '--ploidy', str(args.top_n),
            '--max-projection-gap', str(args.max_projection_gap),
            '--max-projection-blocks', str(args.max_projection_blocks),
            '--min-projection-aligned-bp', str(args.min_projection_aligned_bp),
            '--min-projection-coverage', str(args.min_projection_coverage),
            '--min-block-identity', str(args.min_block_identity),
            '--min-block-mapq', str(args.min_block_mapq),
            '--min-block-collinearity', str(args.min_block_collinearity),
            '--min-anchor-block-aligned-bp', str(args.min_anchor_block_aligned_bp),
            '--min-anchor-query-fraction', str(args.min_anchor_query_fraction),
            '--min-segment-length', str(args.min_segment_length),
            '--min-path-envelope-query-length', str(args.min_path_envelope_query_length),
            '--min-path-envelope-query-coverage', str(args.min_path_envelope_query_coverage),
            '--min-path-envelope-query-span-coverage', str(args.min_path_envelope_query_span_coverage),
            '--min-path-envelope-target-coverage', str(args.min_path_envelope_target_coverage),
            '--max-path-envelope-span-ratio', str(args.max_path_envelope_span_ratio),
            '--min-long-path-pair-overlap', str(args.min_long_path_pair_overlap),
            '--min-over-capacity-pair-overlap', str(args.min_over_capacity_pair_overlap),
            '--min-over-capacity-pair-short-coverage', str(args.min_over_capacity_pair_short_coverage),
            '--min-deferred-over-capacity-bp', str(args.min_deferred_over_capacity_bp),
            '--max-deferred-over-capacity-query-length', str(args.max_deferred_over_capacity_query_length),
            '--max-deferred-over-capacity-preferred-fraction', str(args.max_deferred_over_capacity_preferred_fraction),
            '--length-prior-scale', str(args.length_prior_scale),
            '--unknown-dosage-policy', args.unknown_dosage_policy,
            '--over-capacity-policy', args.over_capacity_policy,
            '--min-resolution-margin', str(args.min_resolution_margin)
        ]
        if args.gfa:
            table_command.extend(['--gfa', args.gfa])
        subprocess.run(table_command, check=True)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

    if stop_after == 'table':
        print(f'PHap v2 allelic table: {allelic_table_file}')
        return

    ### Step 2: Extract chromosome sequences from p_utg
    step2_dir = os.path.join(output_root, '02.chr_seq')
    os.makedirs(step2_dir, exist_ok=True)
    try:
        subprocess.run([
            sys.executable,
            os.path.join(utils_realpath, 'extract_chr_from_putg.py'),
            '--p_utg', args.p_utg,
            '--mT2T', args.mT2T,
            '--paf', best_paf_file,
            '--chain-qc', chain_qc_file,
            '--wd', step2_dir,
            '--chr_num', str(args.chr_num),
            '--min-chromosome-margin', str(args.min_chromosome_margin),
        ], check=True)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

    if stop_after == 'chr_seq':
        print(f'PHap v2 chromosome sequences: {step2_dir}')
        return

    ### step3: Cluster based on hi-c links
    step3_dir = os.path.join(output_root, '03.cluster')
    os.makedirs(step3_dir, exist_ok=True)
    # 待聚类的染色体文件列表
    file_list = sorted(glob.glob(os.path.join(step2_dir, '*.putg.fa')))
    print(f'Debugs: file_list {file_list}')

    for file in file_list:
        file_bn = os.path.basename(file)
        chr = file_bn.replace('.putg.fa', '')
        # 聚类
        cluster_chr_dir = os.path.join(step3_dir, chr)
        os.makedirs(cluster_chr_dir, exist_ok=True)
        os.chdir(cluster_chr_dir)
        allelic_table_chr = f'{chr}.corrected_allelic_table.txt'
        with open(allelic_table_file) as source, open(allelic_table_chr, 'w') as destination:
            for row in source:
                if row.split('\t', 1)[0] == chr:
                    destination.write(row)
        cluster(
            file,
            args.full_links,
            args.flank,
            args.contig_type,
            allelic_table_chr,
            cluster_chr_dir,
            ploidy=args.top_n,
            balance_weight=args.cluster_balance_weight,
            max_refinement_rounds=args.cluster_refinement_rounds,
            max_phase_block_rounds=args.cluster_phase_block_rounds,
            min_phase_block_gain=args.cluster_min_phase_block_gain,
            max_phase_boundary_relaxations=args.cluster_max_phase_boundary_relaxations,
            max_phase_boundary_overlap=args.cluster_max_phase_boundary_overlap,
            max_constraint_relaxation_rounds=args.cluster_max_constraint_relaxation_rounds,
            min_constraint_relaxation_gain=args.cluster_min_constraint_relaxation_gain,
            min_constraint_relaxation_links=args.cluster_min_constraint_relaxation_links,
            min_constraint_relaxation_margin=args.cluster_min_constraint_relaxation_margin,
            max_constraint_relaxations_per_move=args.cluster_max_constraint_relaxations_per_move,
            max_constraint_relaxation_overlap=args.cluster_max_constraint_relaxation_overlap,
            max_phase_interval_rounds=args.cluster_max_phase_interval_rounds,
            max_phase_interval_relaxations=args.cluster_max_phase_interval_relaxations,
            max_backtracks=args.cluster_max_backtracks,
            constraint_relaxation=args.cluster_constraint_relaxation,
            protected_long_unitig_length=args.cluster_protected_long_unitig_length,
            protected_length_ratio=args.cluster_protected_length_ratio,
            protected_short_overlap=args.cluster_protected_short_overlap,
            hic_link_normalization=args.hic_link_normalization,
        )
        # cluster(args.p_utg, args.full_links, args.flank, args.contig_type, allelic_table_chr, cluster_chr_dir)
    os.chdir(cwd)

    if stop_after == 'cluster':
        print(f'PHap v2 cluster results: {step3_dir}')
        return

    ### step4: Re-cluster unclustered chromosome unitigs
    jobs = []
    step4_dir = os.path.join(output_root, '04.recluster')
    os.makedirs(step4_dir, exist_ok=True)
    for file in file_list:
        file_bn = os.path.basename(file)
        chr = file_bn.replace('.putg.fa', '')
        recluster_chr_dir = os.path.join(step4_dir, chr)
        os.makedirs(recluster_chr_dir, exist_ok=True)
        command = [
            sys.executable,
            os.path.join(utils_realpath, 'chr_uncluster_recluster.py'),
            '--fasta', file,
            '--contig-type', args.contig_type,
            '--full-links', args.full_links,
            '--clusters-file', os.path.join(step3_dir, chr, 'group.cluster.txt'),
            '--allelic-table', os.path.join(
                step3_dir, chr, f'{chr}.corrected_allelic_table.txt'
            ),
            '--relaxed-constraints', os.path.join(
                step3_dir, chr, 'cluster_relaxed_constraints.tsv'
            ),
            '--output-dir', recluster_chr_dir,
            '--ploidy', str(args.top_n),
            '--hic-link-normalization', args.hic_link_normalization,
            '--min-adjusted-links', str(args.recluster_min_adjusted_links),
            '--min-group-margin', str(args.recluster_min_group_margin),
            '--min-allelic-block-anchors', str(args.recluster_min_allelic_block_anchors),
            '--max-allelic-block-configurations', str(
                args.recluster_max_allelic_block_configurations
            ),
            '--max-allelic-block-movable-length', str(
                args.recluster_max_allelic_block_movable_length
            ),
            '--min-group-bp-ratio', str(args.recluster_min_group_bp_ratio),
            '--min-assigned-fraction', str(args.recluster_min_assigned_fraction),
            '--max-rounds', str(args.recluster_max_rounds),
            '--refinement-rounds', str(args.recluster_refinement_rounds),
            '--low-confidence-policy', args.recluster_low_confidence_policy,
            '--unknown-dosage-policy', args.recluster_unknown_dosage_policy,
            '--reviewed-seed-fallback', args.recluster_reviewed_seed_fallback,
        ]
        if args.recluster_seed_review == 'weak':
            command.extend([
                '--cluster-assignments',
                os.path.join(step3_dir, chr, 'cluster_assignments.tsv'),
                '--trusted-seed-bases', args.recluster_trusted_seed_bases,
            ])
        if args.flank is not None:
            command.extend(['--flank', str(args.flank)])
        jobs.append((
            command,
            os.path.join(recluster_chr_dir, 'log_re_out'),
            os.path.join(recluster_chr_dir, 'log_re_err'),
        ))
    workers = max(1, min(args.threads, len(jobs)))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(run_logged_command, *job): job[0]
            for job in jobs
        }
        for future in as_completed(futures):
            future.result()

    os.chdir(cwd)

    if stop_after == 'recluster':
        print(f'PHap v2 recluster results: {step4_dir}')
        return

    ### step5: rescue chromosome-unassigned unitigs
    step5_dir = os.path.join(output_root, '05.rescue')
    os.makedirs(step5_dir, exist_ok=True)
    file_un_chr_fasta = os.path.join(step2_dir, 'un_chr.fa')
    rescue_command = [
        sys.executable,
        os.path.join(utils_realpath, 'unchr_recluster.py'),
        '--assembly-fasta', args.p_utg,
        '--candidate-fasta', file_un_chr_fasta,
        '--contig-type', args.contig_type,
        '--full-links', args.full_links,
        '--recluster-dir', step4_dir,
        '--output-dir', step5_dir,
        '--ploidy', str(args.top_n),
        '--hic-link-normalization', args.hic_link_normalization,
        '--min-adjusted-links', str(args.rescue_min_adjusted_links),
        '--min-chromosome-margin', str(args.rescue_min_chromosome_margin),
        '--min-group-margin', str(args.rescue_min_group_margin),
        '--min-assigned-fraction', str(args.rescue_min_assigned_fraction),
        '--max-rounds', str(args.rescue_max_rounds),
        '--low-confidence-policy', args.rescue_low_confidence_policy,
        '--unknown-dosage-policy', args.rescue_unknown_dosage_policy,
    ]
    if args.flank is not None:
        rescue_command.extend(['--flank', str(args.flank)])
    run_logged_command(
        rescue_command,
        os.path.join(step5_dir, 'log_rescue_out'),
        os.path.join(step5_dir, 'log_rescue_err'),
    )

    os.chdir(cwd)


if __name__ == '__main__':
    main()
