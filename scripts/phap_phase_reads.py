#!/usr/bin/env python
'''
time: 2024-07-29
author: pxxiao
version: 1.0
description:
根据 HiFi reads 比对的 bam 文件；contig type；unitigs groups：
    根据 contig 的类型，将 collapsed unitigs HiFi reads 随机分为 2 份、3 份、四份；
    最后得到每个 Group 包含的 HiFi reads；
    然后，根据 HiFi reads，对每个 Group 分别进行 hifiasm 组装。
----
time: 2024-07-30
author: pxxiao
version: 2.0
description:
增加 ONT reads，Hi-C reads 的 group-specific reads 划分；
然后，使用 HiFi+ONT 策略进行组装；
haphic 利用 Hi-C reads 进行挂载（还得看一下效果如何）
'''


import argparse
import sys
from pathlib import Path

from phap_core.runner import PreflightError, require_input_files, run_command


def parse_args():
    parser = argparse.ArgumentParser("Haplotype assembly and scaffolding of autopolyploid genome.")
    parser.add_argument('--bam_hifi', type=str, required=True, help='Path to HiFi bam file')
    parser.add_argument('--bam_hic', type=str, required=True, help='Path to Hi-C bam file')
    parser.add_argument('--bam_ont', type=str, required=True, help='Path to ONT bam file')
    parser.add_argument('--contig_type', required=True, type=str, help='Path to contig type from dosage analysis')
    parser.add_argument('--group', type=str, help='Path to merge group file from cluster, default is /02.cluster/05.rescue/merge.group.reassignment.cluster.txt')
    parser.add_argument('--hifi', type=str, required=True, help='Path to HiFi reads file')
    parser.add_argument('--ont', type=str, required=True, help='Path to ONT reads file')
    parser.add_argument('--hic1', type=str, required=True, help='Path to Hi-C forward reads file')
    parser.add_argument('--hic2', type=str, required=True, help='Path to Hi-C reverse reads file')
    parser.add_argument('--ont_length', type=int, default=1, help='Length of ONT reads file for hifiasm asm [1]')
    parser.add_argument('--ont_quality', type=int, default=0, help='Quality of ONT [0]')
    parser.add_argument('--threads', type=int, default=10, help='The number of threads [10]')
    parser.add_argument('--process', type=int, default=4, help='The number of processes [4]')
    parser.add_argument('--seed', type=int, default=100, help='Random seed for reproducibility [100]')
    parser.add_argument('--ploidy', type=int, required=True, help='Genome ploidy')
    parser.add_argument(
        '--min_mapq',
        type=int,
        default=1,
        help='Minimum MAPQ for read assignment [1]',
    )

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    cwd = Path.cwd()
    group_file = (
        Path(args.group).expanduser().resolve()
        if args.group
        else cwd / '02.cluster' / '05.rescue' / 'group.reassignment.cluster.txt'
    )
    input_paths = [
        args.bam_hifi,
        args.bam_hic,
        args.bam_ont,
        args.contig_type,
        group_file,
        args.hifi,
        args.ont,
        args.hic1,
        args.hic2,
    ]
    try:
        require_input_files(input_paths)
    except PreflightError as error:
        raise SystemExit(f"phap phase_reads: error: {error}") from error

    ### Step 1: p_utg vs mT2T
    step1_dir = cwd / '03.phase_reads'
    step1_dir.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        '-m',
        'utils.phase_reads_assemble_anchor',
        '--bam_hifi',
        str(Path(args.bam_hifi).expanduser().resolve()),
        '--bam_hic',
        str(Path(args.bam_hic).expanduser().resolve()),
        '--bam_ont',
        str(Path(args.bam_ont).expanduser().resolve()),
        '--contig_type',
        str(Path(args.contig_type).expanduser().resolve()),
        '--group',
        str(group_file),
        '--hifi',
        str(Path(args.hifi).expanduser().resolve()),
        '--ont',
        str(Path(args.ont).expanduser().resolve()),
        '--ont_length',
        str(args.ont_length),
        '--ont_quality',
        str(args.ont_quality),
        '--hic1',
        str(Path(args.hic1).expanduser().resolve()),
        '--hic2',
        str(Path(args.hic2).expanduser().resolve()),
        '--threads',
        str(args.threads),
        '--process',
        str(args.process),
        '--seed',
        str(args.seed),
        '--ploidy',
        str(args.ploidy),
        '--min_mapq',
        str(args.min_mapq),
    ]
    with (step1_dir / 'log_phase_reads_assemble_anchor_out').open('w') as stdout, (
        step1_dir / 'log_phase_reads_assemble_anchor_err'
    ).open('w') as stderr:
        run_command(command, cwd=step1_dir, stdout=stdout, stderr=stderr)


if __name__ == '__main__':
    main()
