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
----
time: 2024-09-29
author: pxxiao
version: 3.0
description:
添加一个参数，控制 ONT reads 的长度，用于 hifiasm UL 组装
'''


import argparse
import logging
import os
import shlex
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Tuple

import pysam  # type: ignore[import-not-found]

from phap_core.read_assignment import (
    ReadAssignment,
    UnitigCandidate,
    alignment_is_usable,
    assign_reads,
    build_unitig_candidates,
    canonical_read_id,
    group_assigned_reads,
    paired_fastq_patterns,
    parse_unitig_dosages,
)
from phap_core.runner import (
    PreflightError,
    require_input_files,
    require_tools,
    run_shell_command,
    run_shell_commands_parallel,
)

run_in_parallel = run_shell_commands_parallel

# 设置日志记录配置
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


def parse_contig_type_based_on_dosage(
    file_contig_type: str,
) -> Tuple[Dict[str, Optional[int]], Dict[str, str]]:
    """Parse numeric dosage or legacy contig types without guessing invalid calls."""
    with open(file_contig_type) as source:
        dosage_by_unitig, source_labels = parse_unitig_dosages(source)
    logging.info('Parsed dosage for %d unitigs', len(dosage_by_unitig))
    return dosage_by_unitig, source_labels


def parse_group_cluster(group_file: str) -> Dict[str, list[str]]:
    dic_group_ctg: Dict[str, list[str]] = {}
    with open(group_file) as source:
        for lines in source:
            if not lines.strip():
                continue
            line = lines.strip().split()
            group_name = line[0]
            if group_name in dic_group_ctg:
                raise ValueError(f'duplicate group name: {group_name}')
            dic_group_ctg[group_name] = line[1:]
    logging.info('Parsed group clusters successfully')
    return dic_group_ctg


def parse_bam_read_unitigs(
    bam_file: str,
    min_mapq: int,
    paired: bool = False,
) -> Dict[str, set[str]]:
    """Collect one entity per read/pair under the documented BAM filter policy.

    Unmapped, secondary, supplementary, duplicate, QC-fail, and alignments
    below ``min_mapq`` are filtered. Non-proper alignments are intentionally
    retained because that flag is not meaningful evidence against a Hi-C pair.
    """
    read_unitigs: Dict[str, set[str]] = defaultdict(set)
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read.query_name is None:
                raise ValueError(f'alignment without query name in {bam_file}')
            read_id = canonical_read_id(read.query_name, paired=paired)
            read_unitigs[read_id]  # retain filtered/unmapped entities for audit
            if not alignment_is_usable(
                is_unmapped=read.is_unmapped,
                is_secondary=read.is_secondary,
                is_supplementary=read.is_supplementary,
                is_duplicate=read.is_duplicate,
                is_qcfail=read.is_qcfail,
                is_proper_pair=read.is_proper_pair,
                mapping_quality=read.mapping_quality,
                min_mapq=min_mapq,
            ):
                continue
            contig = bam.get_reference_name(read.reference_id)
            if contig is not None:
                read_unitigs[read_id].add(contig)
    logging.info('Parsed %d read entities from %s', len(read_unitigs), bam_file)
    return read_unitigs


def atomic_write_lines(path: str, lines: Iterable[str]) -> None:
    target = Path(path)
    temporary = target.with_name(target.name + '.tmp')
    with temporary.open('w') as output:
        for line in lines:
            output.write(line)
            output.write('\n')
    os.replace(temporary, target)


def write_assignment_audits(
    candidates: Mapping[str, UnitigCandidate],
    source_labels: Mapping[str, str],
    assignments: Iterable[ReadAssignment],
    groups: Iterable[str],
    *,
    seed: int,
    ploidy: int,
    min_mapq: int,
    filter_policy: str,
) -> None:
    assignment_list = sorted(assignments, key=lambda item: (item.modality, item.read_id))
    candidate_lines = [
        'unitig_ID\tdosage\tcandidate_groups\tstatus\treason\tsource_label'
    ]
    for unitig in sorted(candidates):
        candidate = candidates[unitig]
        dosage = '' if candidate.dosage is None else str(candidate.dosage)
        candidate_lines.append(
            '\t'.join(
                [
                    candidate.unitig_id,
                    dosage,
                    ','.join(candidate.groups),
                    candidate.status,
                    candidate.reason,
                    source_labels.get(unitig, 'missing'),
                ]
            )
        )
    atomic_write_lines('unitig_candidates.tsv', candidate_lines)

    read_lines = [
        'read_ID\tmodality\tstatus\tdestination_group\tcandidate_groups\tunitigs\treason'
    ]
    for assignment in assignment_list:
        read_lines.append(
            '\t'.join(
                [
                    assignment.read_id,
                    assignment.modality,
                    assignment.status,
                    assignment.destination_group or '',
                    ','.join(assignment.candidate_groups),
                    ','.join(assignment.unitigs),
                    assignment.reason,
                ]
            )
        )
    atomic_write_lines('read_assignments.tsv', read_lines)

    group_names = sorted(set(groups))
    summary_lines = [
        'modality\tgroup\tassigned_reads\tambiguous_reads\tunassigned_reads\t'
        'total_entities\tcross_group_overlaps\tcross_group_overlap_rate\t'
        'estimated_cross_group_contamination_rate\tseed\tploidy\tmin_mapq\t'
        'filter_policy'
    ]
    modalities = sorted({assignment.modality for assignment in assignment_list})
    for modality in modalities:
        subset = [item for item in assignment_list if item.modality == modality]
        ambiguous = sum(item.status == 'ambiguous' for item in subset)
        unassigned = sum(item.status == 'unassigned' for item in subset)
        assigned_ids = [
            item.read_id for item in subset if item.status == 'assigned'
        ]
        owner_counts: Dict[str, int] = defaultdict(int)
        for read_id in assigned_ids:
            owner_counts[read_id] += 1
        overlaps = sum(count > 1 for count in owner_counts.values())
        overlap_rate = overlaps / len(assigned_ids) if assigned_ids else 0.0
        for group in group_names:
            assigned = sum(
                item.status == 'assigned' and item.destination_group == group
                for item in subset
            )
            summary_lines.append(
                f'{modality}\t{group}\t{assigned}\t0\t0\t{assigned}\t0\t0.000000\t'
                f'0.000000\t{seed}\t{ploidy}\t{min_mapq}\t{filter_policy}'
            )
        summary_lines.append(
            f'{modality}\tALL\t{len(assigned_ids)}\t{ambiguous}\t{unassigned}\t'
            f'{len(subset)}\t{overlaps}\t{overlap_rate:.6f}\t'
            f'{overlap_rate:.6f}\t{seed}\t{ploidy}\t{min_mapq}\t{filter_policy}'
        )
    atomic_write_lines('read_assignment_summary.tsv', summary_lines)


def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(description="Haplotype assembly and scaffolding of autopolyploid genome")
    parser.add_argument('--bam_hifi', type=str, required=True, help='Path to HiFi bam file')
    parser.add_argument('--bam_hic', type=str, required=True, help='Path to Hi-C bam file')
    parser.add_argument('--bam_ont', type=str, required=True, help='Path to ONT bam file')
    parser.add_argument('--contig_type', required=True, type=str, help='Path to contig type from dosage analysis')
    parser.add_argument('--group', type=str, required=True, help='Path to group file from cluster')
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
        help='Minimum MAPQ for primary, mapped, non-duplicate alignments [1]',
    )

    args = parser.parse_args()
    return args


def main() -> None:
    args = parse_args()
    if args.ploidy < 1:
        raise SystemExit('phap phase_reads: error: --ploidy must be at least 1')
    if args.min_mapq < 0:
        raise SystemExit('phap phase_reads: error: --min_mapq must be non-negative')
    required_tools = [
        'awk',
        'bash',
        'bwa',
        'filter_bam',
        'haphic',
        'hifiasm',
        'pigz',
        'samblaster',
        'samtools',
        'seqkit',
    ]
    if args.ont_length != 1 or args.ont_quality != 0:
        required_tools.append('chopper')
    try:
        require_input_files([
            args.bam_hifi,
            args.bam_hic,
            args.bam_ont,
            args.contig_type,
            args.group,
            args.hifi,
            args.ont,
            args.hic1,
            args.hic2,
        ])
        require_tools(required_tools)
    except PreflightError as error:
        raise SystemExit(f"phap phase_reads: error: {error}") from error

    try:
        ### steps 1-3: parse inputs and make one global decision per read entity
        dic_contig_dosage, source_labels = parse_contig_type_based_on_dosage(
            args.contig_type
        )
        dic_group_ctg = parse_group_cluster(args.group)
        candidates = build_unitig_candidates(
            dic_group_ctg,
            dic_contig_dosage,
            ploidy=args.ploidy,
        )
        bam_inputs = {
            'hifi': (args.bam_hifi, False),
            'ont': (args.bam_ont, False),
            'hic': (args.bam_hic, True),
        }
        assignments_by_modality: Dict[str, Tuple[ReadAssignment, ...]] = {}
        all_assignments: list[ReadAssignment] = []
        for modality in sorted(bam_inputs):
            bam_file, paired = bam_inputs[modality]
            read_unitigs = parse_bam_read_unitigs(
                bam_file,
                min_mapq=args.min_mapq,
                paired=paired,
            )
            assignments = assign_reads(
                read_unitigs,
                candidates,
                modality=modality,
                seed=args.seed,
            )
            assignments_by_modality[modality] = assignments
            all_assignments.extend(assignments)

        groups = sorted(dic_group_ctg)
        dic_group_hifi = group_assigned_reads(assignments_by_modality['hifi'], groups)
        dic_group_ont = group_assigned_reads(assignments_by_modality['ont'], groups)
        dic_group_hic = group_assigned_reads(assignments_by_modality['hic'], groups)
        filter_policy = (
            'filter=unmapped,secondary,supplementary,duplicate,qcfail,'
            'mapq_below_min;retain=primary_mapped_including_nonproper'
        )
        write_assignment_audits(
            candidates,
            source_labels,
            all_assignments,
            groups,
            seed=args.seed,
            ploidy=args.ploidy,
            min_mapq=args.min_mapq,
            filter_policy=filter_policy,
        )
    except (OSError, ValueError) as error:
        raise SystemExit(f'phap phase_reads: error: {error}') from error

    ont_reads = args.ont
    ont_reads_bn = os.path.basename(ont_reads)
    ont_quality = args.ont_quality
    ont_length = args.ont_length
    if args.ont_length != 1 or args.ont_quality != 0:
        ont_reads_filter = ont_reads_bn + '.filter.ont.fq'
        ont_reads_filter_gzip = ont_reads_bn + '.filter.ont.fq.gz'
        cmd = (
            f"chopper -l {ont_length} -q {ont_quality} -t {args.threads} "
            f"-i {shlex.quote(ont_reads)} > {shlex.quote(ont_reads_filter)} && "
            f"pigz -p {args.threads} {shlex.quote(ont_reads_filter)}"
        )
        run_shell_command(cmd)

    commands = []
    process_seqkit = args.process
    for group, reads in dic_group_hifi.items():
        atomic_write_lines(f'{group}.HiFi.txt', reads)
        cmd = f'seqkit grep -j {args.threads} -f {group}.HiFi.txt {args.hifi} > {group}.HiFi.fq && pigz -p {args.threads} {group}.HiFi.fq'
        commands.append(cmd)
        logging.info(f'Executing command: {cmd}')
    for group, reads in dic_group_ont.items():
        atomic_write_lines(f'{group}.ONT.txt', reads)
        if args.ont_length != 1 or args.ont_quality != 0:
            cmd = f'seqkit grep -j {args.threads} -f {group}.ONT.txt {ont_reads_filter_gzip} > {group}.ONT.fq && pigz -p {args.threads} {group}.ONT.fq'
        else:
            cmd = f'seqkit grep -j {args.threads} -f {group}.ONT.txt {args.ont} > {group}.ONT.fq && pigz -p {args.threads} {group}.ONT.fq'
        commands.append(cmd)
        logging.info(f'Executing command: {cmd}')
    for group, reads in dic_group_hic.items():
        mate1_patterns = sorted(
            {
                pattern
                for read in reads
                for pattern in paired_fastq_patterns(read, mate=1)
            }
        )
        mate2_patterns = sorted(
            {
                pattern
                for read in reads
                for pattern in paired_fastq_patterns(read, mate=2)
            }
        )
        atomic_write_lines(f'{group}.Hi-C.1.txt', mate1_patterns)
        atomic_write_lines(f'{group}.Hi-C.2.txt', mate2_patterns)
        cmd1 = f'seqkit grep -j {args.threads} -f {group}.Hi-C.1.txt {args.hic1} > {group}.Hi-C.1.fq && pigz -p 5 {group}.Hi-C.1.fq'
        cmd2 = f'seqkit grep -j {args.threads} -f {group}.Hi-C.2.txt {args.hic2} > {group}.Hi-C.2.fq && pigz -p 5 {group}.Hi-C.2.fq'
        commands.append(cmd1)
        commands.append(cmd2)
        logging.info(f'Executing command: {cmd1}')
        logging.info(f'Executing command: {cmd2}')
    run_in_parallel(commands, process_seqkit)

    ### step5: 对每个 group 进行组装（hifiasm）
    commands = []
    for group in dic_group_hifi.keys():
        os.makedirs(f'{group}.asm', exist_ok=True)
        cmd_hifiasm = f'hifiasm -t {args.threads} -o {group}.asm/{group}.asm --ul {group}.ONT.fq.gz {group}.HiFi.fq.gz > {group}.asm/log_hifiasm_out 2> {group}.asm/log_hifiasm_err'
        commands.append(cmd_hifiasm)
        logging.info(f'Executing command: {cmd_hifiasm}')
    run_in_parallel(commands, args.process)

    ### step6: 对每个 group p_ctg 进行挂载（haphic）
    commands = []
    for group in dic_group_hifi.keys():
        # 创建目录
        os.makedirs(f'{group}.asm/scaffolding/01_hic_mapping', exist_ok=True)
        os.makedirs(f'{group}.asm/scaffolding/02_haphic', exist_ok=True)
        # 构建命令序列
        cmd_scaffolding = [
            f'awk \'{{if($0~/^S/) print ">"$2"\\n"$3}}\' {group}.asm/{group}.asm.bp.p_ctg.gfa > {group}.asm/scaffolding/01_hic_mapping/{group}.asm.bp.p_ctg.gfa.fa',
            f'bwa index {group}.asm/scaffolding/01_hic_mapping/{group}.asm.bp.p_ctg.gfa.fa',
            f'bwa mem -5SP -t {args.threads} {group}.asm/scaffolding/01_hic_mapping/{group}.asm.bp.p_ctg.gfa.fa {group}.Hi-C.1.fq.gz {group}.Hi-C.2.fq.gz | samblaster | samtools view - -@ {args.threads} -S -h -b -F 3340 -o {group}.asm/scaffolding/01_hic_mapping/HiC.bam',
            f'filter_bam {group}.asm/scaffolding/01_hic_mapping/HiC.bam 1 --nm 3 --threads {args.threads} | samtools view - -b -@ {args.threads} -o {group}.asm/scaffolding/01_hic_mapping/HiC.filtered.bam',
            f'cd {group}.asm/scaffolding/02_haphic',
            f'haphic pipeline ../01_hic_mapping/{group}.asm.bp.p_ctg.gfa.fa ../01_hic_mapping/HiC.filtered.bam 1 --threads {args.threads} --processes {args.process} --Nx 100 > log_haphic_out 2> log_haphic_err',
            'cd 04.build',
            'bash juicebox.sh > log_juicebox_out 2> log_juicebox_err'
        ]
        # 将命令加入到 commands 列表中
        commands.append(" && ".join(cmd_scaffolding))
        logging.info(f'Executing command sequence for {cmd_scaffolding}')
    run_in_parallel(commands, args.process)


if __name__ == '__main__':
    main()
