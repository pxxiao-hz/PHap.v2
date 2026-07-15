#!/usr/bin/env python

import os
from collections import defaultdict

from phap_core.runner import run_shell_command, run_shell_commands_parallel

run_in_parallel = run_shell_commands_parallel

### 制作文件夹，先判断文件夹是否存在？
def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    return None


def check_file_in_path(file, cmd):
    if os.path.exists(file):
        print(f'[info] {file} exists, CMD: {cmd}; PASS!')
    else:
        run_shell_command(cmd)
    return None


### 判断文件是否存在，然后执行某个函数
def check_file_in_path_fun(file, fun):
    '''
    文件不存在，执行某个函数；存在则跳过
    :param file:
    :param fun:
    :return:
    '''
    if os.path.exists(file):
        print('File exits')
    else:
        print('File not exits')
        # fun()
    return None


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


def delta_filter_same_seq(file_delta, file_filtered_delta):
    '''
    对于short to short比对，需要去除比对到自身的部分
    :param file_delta: 待过滤的 delta 文件
    :param file_filtered_delta: 过滤相同 seq 匹配之后的 delta 文件
    :return:
    '''
    out = open(file_filtered_delta, 'w')
    for lines in open(file_delta, 'r'):
        line = lines.strip().split()
        if lines.startswith('/'):
            out.write(lines)
            continue
        if lines.startswith('NUCMER'):
            out.write(lines)
            continue
        ## >ptg000011l ptg000028l 41964882 2810910
        if lines.startswith('>'):
            if line[0][1:] != line[1]:  # 不是自身比对
                block = 1
                out.write(lines)
            else:
                block = 0               # 自身比对
        ## 具体的比对信息
        if len(line) == 1 or len(line) == 7:
            if block == 0:
                continue
            else:
                out.write(lines)
    return None


def cal_aln_rate(paf_file):
    def interval_merge(list):
        # list: [a,b], [c,d], [e,f]...
        merged = []
        # 对list进行排序
        sorted_list = sorted(list, key=lambda x: x[0])
        for interval in sorted_list:
            if not merged or interval[0] > merged[-1][1]:  # 判断是不是第一个区间，或者当前区间的开始是否大于已有区间的结尾
                merged.append(interval)
            else:  # 当前区间的开始小于已有区间的结尾，进行合并
                merged[-1][1] = max(merged[-1][1], interval[1])
        return merged

    def sum_interval(list):
        su = 0
        for i in list:
            su += (i[1] - i[0])
        return su

    output_file = paf_file + 'aln.ratio.txt'
    output_file = open(output_file, 'w')
    dic_paf = parse_paf(paf_file)
    print('qryID\trefID\tqryLen\trefLen\taln_direction\tratio_qry\tratio_ref\tratio_qry2ref')
    output_file.write('#qryID\trefID\tqryLen\trefLen\taln_direction\tratio_qry\tratio_ref\tratio_qry2ref\n')
    for key, values in dic_paf.items():
        qry_len = values[0][0]
        ref_len = values[1][0]
        # print(qry_len, ref_len)
        ref_list = interval_merge(values[3])
        qry_list = interval_merge(values[2])
        positive_qry_list = interval_merge(values[4])
        negative_qry_list = interval_merge(values[5])
        sum_ref = sum_interval(ref_list)  # 比对上的长度：ref
        sum_qry = sum_interval(qry_list)  # 比对上的长度：qry
        sum_positive_qry = sum_interval(positive_qry_list)  # 比对上的，正向的长度
        sum_negative_qry = sum_interval(negative_qry_list)  # 比对上的，负向的长度

        # 各自比对上的比例
        ratio_ref = sum_ref / ref_len
        ratio_qry = sum_qry / qry_len
        # qry / ref 的比例
        ratio_qry2ref = qry_len / ref_len
        # 比对的方向，positive or negative
        if sum_positive_qry > sum_negative_qry:
            aln_direction = '+'
        else:
            aln_direction = '-'

        ref_ID = key.split('_to_')[1]
        qry_ID = key.split('_to_')[0]
        print(f'{qry_ID}\t{ref_ID}\t{qry_len}\t{ref_len}\t{aln_direction}\t{ratio_qry}\t{ratio_ref}\t{ratio_qry2ref}')
        output_file.write(f'{qry_ID}\t{ref_ID}\t{qry_len}\t{ref_len}\t{aln_direction}\t{ratio_qry}\t{ratio_ref}\t{ratio_qry2ref}\n')
    return dic_paf


def parse_paf(paf_file):
    ''' 解析 paf 文件 '''
    # 字典：存储长度信息
    dic_paf = {}
    for lines in open(paf_file, "r"):
        line = lines.strip().split()
        qry_ID, qry_len, qry_start, qry_end, qry_strand = line[0], int(line[1]), int(line[2]), int(line[3]), line[4]
        tag_ID, tag_len, tag_start, tag_end = line[5], int(line[6]), int(line[7]), int(line[8])
        match_count, block_length, mapping_quality = line[9], line[10], line[11]

        qry_to_tag = qry_ID + "_to_" + tag_ID
        tag_to_qry = tag_ID + "_to_" + qry_ID
        if tag_to_qry not in dic_paf:       # 20240416 fix bugs: ctg1_VS_ctg2 与 ctg2_VS_ctg1只保留一个比对
            if qry_to_tag not in dic_paf.keys():
                # qry len, tag len, qry coordinate interval, tag coordinate interval, sum_qry_forward, sum_qry_reverse
                dic_paf[qry_to_tag] = [[], [], [], [], [], []]
                dic_paf[qry_to_tag][0].append(qry_len)
                dic_paf[qry_to_tag][1].append(tag_len)
                dic_paf[qry_to_tag][2].append([qry_start, qry_end])
                dic_paf[qry_to_tag][3].append([tag_start, tag_end])
            else:
                dic_paf[qry_to_tag][2].append([qry_start, qry_end])
                dic_paf[qry_to_tag][3].append([tag_start, tag_end])
            if qry_strand == '+':
                # 链一致
                dic_paf[qry_to_tag][4].append([qry_start, qry_end])
            elif qry_strand == '-':
                # 链不一致
                dic_paf[qry_to_tag][5].append([qry_start, qry_end])
    print('20240416 fix bugs: 过滤之前，互相比对的染色体对数', len(dic_paf))
    return dic_paf


def parse_delta(file_delta):

    ### delta's basename
    bn_delta = os.path.basename(file_delta).replace('.delta', '')

    ### step1: 过滤 DELTA 文件
    file_filtered_delta = bn_delta + '.filtered.same.seq.delta'
    print('pctg_self_align.filtered.same.seq.delta')
    # delta_filter_same_seq(file_delta, file_filtered_delta)
    check_file_in_path_fun(file_filtered_delta, lambda: delta_filter_same_seq(file_delta, file_filtered_delta))

    ### step2: 将过滤之后的 DELTA 文件转为 PAF 文件
    file_filtered_delta_paf = bn_delta + '.filtered.same.seq.delta.paf'
    cmd_delta2paf = f'delta2paf {file_filtered_delta} > {file_filtered_delta_paf}'
    check_file_in_path(file_filtered_delta_paf, cmd_delta2paf)

    ### step3: 计算比对率
    cal_aln_rate(f'{file_filtered_delta_paf}')

    ### step4: 去除冗余
    list_remove = []
    # for lines in open('pctg_self_align.filtered.same.seq.delta.pafaln.ratio.txt', 'r'):
    for lines in open('pctg_self_align_maxmatch.filtered.same.seq.delta.pafaln.ratio.txt', 'r'):
        if lines.startswith('#'):
            continue
        line = lines.strip().split()
        len_qry = int(line[2])
        len_ref = int(line[3])
        ratio_qry = float(line[5])
        ratio_ref = float(line[6])

        ## 判断 reference 与 query 的长短
        if len_qry < len_ref:
            if ratio_qry > 0.1 and line[0] not in list_remove:
                print(f'Ref length long: {line[0]}')
                print(line)
                list_remove.append(line[0])
        elif len_ref < len_qry:
            if ratio_ref > 0.1 and line[1] not in list_remove:
                print(f'Qry length long: {line[1]}')
                print(line)
                list_remove.append(line[1])

    for i in list_remove:
        print(i)
    return None


def process_delta(file_delta, args):

    def delta2paf(file_delta, file_paf):
        ''' convert delta to PAF '''
        cmd_delta2paf = f'delta2paf {file_delta} > {file_paf}'
        check_file_in_path(file_paf, cmd_delta2paf)

    def get_contig_pair(contig1, contig2, len_contig1, len_contig2):
        ''' contig pair based on contig length '''
        if len_contig1 > len_contig2:
            return (contig1, contig2)
        else:
            return (contig2, contig1)

    def filter_paf(file_paf, file_paf_filter, length_alignment, length_query):
        ''' filter PAF '''
        dic_paf_info = defaultdict(lambda: [0, 0, 0, 0, 0.00, 0.00])
        processed_contig_pair = set()
        for lines in open(file_paf, 'r'):
            line = lines.strip().split()
            ## same seq
            if line[0] == line[5]: continue
            ## alignment length
            if int(line[9]) < length_alignment: continue
            ## query length
            if int(line[1]) < length_query: continue
            ## 判断是否处理过？ctg1, ctg2
            if (line[5], line[0]) not in processed_contig_pair:
                processed_contig_pair.add((line[0], line[5]))
            else:
                continue
            ## contig pair
            contig_pair = get_contig_pair(line[0], line[5], int(line[1]), int(line[6]))
            dic_paf_info[contig_pair][0] = int(line[1])     # # ctg1_len, ctg2_len, +=alignment_qry_len, +=alignment_ref_len
            dic_paf_info[contig_pair][1] = int(line[6])
            dic_paf_info[contig_pair][2] += int(line[9])
            dic_paf_info[contig_pair][3] += int(line[10])
            dic_paf_info[contig_pair][4] = dic_paf_info[contig_pair][2] / dic_paf_info[contig_pair][0]
            dic_paf_info[contig_pair][5] = dic_paf_info[contig_pair][3] / dic_paf_info[contig_pair][1]
        return dic_paf_info

    ## delta2paf
    basename = os.path.basename(file_delta).replace('.delta', '')
    file_paf = basename + '.paf'
    delta2paf(file_delta, file_paf)

    ## filter paf
    file_paf_filter = basename + '.filter.paf'
    dic_paf_info = filter_paf(file_paf, file_paf_filter, args.min_alginment_length, args.min_contig_length)

    ## Remove redundancy
    remove_contig = set()
    out = open('removed_contigs_based_on_mummer.id.txt', 'w')
    for contg_pair, values in dic_paf_info.items():
        if values[4] > args.match_ratio or values[5] > args.match_ratio:
            contig_long = contg_pair[0]     # 长的 contig
            contig_short = contg_pair[1]    # 短的 contig
            if contig_short not in remove_contig:
                remove_contig.add(contig_short)
            # print(contg_pair, values)
            # print(contig_short)
    for i in remove_contig:
        out.write(i + '\n')
    return dic_paf_info
