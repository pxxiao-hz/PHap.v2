#!/usr/bin/env python
'''
time: 2024-04-27
author: pxxiao
version: 1.0
description: 根据 reads 深度和可配置的单拷贝深度，将 contig 划分为整数剂量类型

输入文件：cnv_winsize10000_step10000_hq.txt
utg000001l      1       10000   1       6       1.6084  1301636
utg000001l      10001   20000   1       45      13.1233 1301636
utg000001l      20001   30000   1       106     30.727  1301636
utg000001l      30001   40000   1       87      25.3734 1301636
utg000001l      40001   50000   1       63      18.182  1301636
utg000001l      50001   60000   1       67      19.3853 1301636
----
time: 2024-06-24
author: pxxiao
version: 2.0
description: 在浮点数的比较那里出现了问题，将浮点数四舍五入，保留两位小数。
'''


import argparse
import pandas as pd


NAMED_DOSAGE_TYPES = {
    1: "haplotig",
    2: "diplotig",
    3: "triplotig",
    4: "tetraplotig",
    5: "pentaplotig",
    6: "hexaplotig",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Classify contigs by integer copy-number dosage"
    )
    parser.add_argument('--input_file', required=True)
    parser.add_argument('--pandepth', action='store_true', help='input from pandepth or not')
    parser.add_argument(
        '--ploidy', type=int, default=4,
        help='Expected haplotype count; depths above this dosage are repeats [4]'
    )
    parser.add_argument(
        '--base-depth', type=float, default=28.0,
        help='Estimated single-copy (haplotig) depth [28]'
    )
    args = parser.parse_args()
    if args.ploidy < 2:
        parser.error('--ploidy must be at least two')
    if args.base_depth <= 0:
        parser.error('--base-depth must be positive')
    return args


def dosage_type_name(dosage):
    return NAMED_DOSAGE_TYPES.get(dosage, f'dosage_{dosage}')


def classify_contig_normal(average_depth, base_value=28, ploidy=4):
    """Classify depth into the nearest positive integer dosage up to ploidy."""
    average_depth = round(average_depth, 2)
    if average_depth < 0.5 * base_value:
        return 'other'
    dosage = int(average_depth / base_value + 0.5)
    if dosage > ploidy:
        return 'replotig'
    return dosage_type_name(max(dosage, 1))


def main():
    args = parse_args()

    ## 获得每个 contig 的平均深度
    # 读取数据文件
    data = pd.read_csv(args.input_file, sep='\t', header=None)

    # 按照 contig ID 分组计算平均深度
    if args.pandepth:
        contig_depth = data.groupby(0)[7].mean().reset_index()
    else:
        contig_depth = data.groupby(0)[5].mean().reset_index()

    # 为结果添加列名
    contig_depth.columns = ['contig_ID', 'average_depth']

    # 添加 contig 类型列
    contig_depth['contig_type'] = contig_depth['average_depth'].apply(
        lambda depth: classify_contig_normal(depth, args.base_depth, args.ploidy)
    )

    # 打印结果
    print(contig_depth)

    # 将 DataFrame 以文件的形式输出
    contig_depth.to_csv("contig_depth.csv", index=False)
    contig_depth.to_csv("contig_depth.txt", sep="\t", index=False)

    #

    return None


if __name__ == '__main__':
    main()
