#!/usr/bin/env python
"""Build a conflict-aware mosaic T2T reference from primary contigs."""


import argparse
import os
import subprocess
from pathlib import Path

from phap_core.fasta import read_fasta
from phap_core.mt2t_overlap import Mt2tOverlapParameters
from phap_core.runner import PreflightError, require_input_files, require_tools, run_command
from utils.mt2t_alignment_workflow import build_pairwise_overlap_paf
from utils.mt2t_overlap_workflow import run_mt2t_overlap_workflow


def parse_args():
    parser = argparse.ArgumentParser('Get mosaic T2T (mT2T) reference from primary contig assembly (p_ctg).')
    # Parameters for mT2T
    mT2T = parser.add_argument_group('>>> mT2T reference')
    mT2T.add_argument('--p_ctg', required=True,
                                        help='path to p_ctg file')
    mT2T.add_argument('--min_distance', required=False, type=float, default=0.15,
                                        help='minimum distance between two contigs [0.15]')
    mT2T.add_argument('--min-longer-coverage', '--match_ratio_target',
                      dest='min_longer_coverage', type=float, default=0.05,
                      help='minimum interval-union coverage of the longer contig [0.05]')
    mT2T.add_argument('--min-shorter-coverage', '--match_ratio_query',
                      dest='min_shorter_coverage', type=float, default=0.1,
                      help='minimum interval-union coverage of the shorter contig [0.1]')
    mT2T.add_argument('--min-identity', type=float, default=0.8,
                      help='minimum alignment-chain identity [0.8]')
    mT2T.add_argument('--min_contig_length', required=False, type=int, default=100000,
                                        help='minimum contig length [100000]')
    mT2T.add_argument('--min_alignment_length', required=False, type=int, default=200,
                                        help='minimum alignment length [200]')
    parser.add_argument('--internal_margin_ratio', type=float, default=0.05,
                        help='Proportional margin to define internal alignments (default: 0.05)')
    parser.add_argument('--max-query-gap', type=int, default=3000000,
                        help='Maximum query gap in a validated overlap chain [3000000]')
    parser.add_argument('--max-target-gap', type=int, default=3000000,
                        help='Maximum oriented target gap in a validated overlap chain [3000000]')
    parser.add_argument('--output-directory', default='01.mT2T',
                        help='Output directory [01.mT2T]')

    # global
    parser.add_argument('--threads', required=False, type=int, default=1, help='number of threads [1]')
    parser.add_argument('--process', required=False, type=int, default=1,
                        help='maximum concurrent external processes [1]')
    parser.add_argument('--cpu-budget', type=int, default=os.cpu_count() or 1,
                        help='total CPU budget across concurrent alignments')

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    try:
        if args.threads < 1 or args.process < 1 or args.cpu_budget < 1:
            raise ValueError('--threads, --process, and --cpu-budget must be positive')
        if args.threads * args.process > args.cpu_budget:
            raise ValueError(
                '--threads * --process exceeds --cpu-budget; reduce concurrency '
                'or explicitly increase the budget'
            )
        if not 0.0 <= args.min_distance <= 1.0:
            raise ValueError('--min_distance must be in [0, 1]')
        parameters = Mt2tOverlapParameters(
            min_contig_length=args.min_contig_length,
            min_alignment_block_length=args.min_alignment_length,
            min_identity=args.min_identity,
            min_shorter_coverage=args.min_shorter_coverage,
            min_longer_coverage=args.min_longer_coverage,
            internal_margin_ratio=args.internal_margin_ratio,
            max_query_gap=args.max_query_gap,
            max_target_gap=args.max_target_gap,
        )
        require_input_files([args.p_ctg])
        fasta_records = read_fasta(args.p_ctg)
        eligible_count = sum(
            len(record.sequence) >= args.min_contig_length
            for record in fasta_records
        )
        if eligible_count >= 2:
            require_tools(["mash", "minimap2"])
    except (PreflightError, ValueError) as error:
        raise SystemExit(f"phap mt2t: error: {error}") from error

    cwd = Path.cwd()
    output_root = Path(args.output_directory)
    if not output_root.is_absolute():
        output_root = cwd / output_root

    tool_versions = (
        {
            'mash': _tool_version('mash'),
            'minimap2': _tool_version('minimap2'),
        }
        if eligible_count >= 2
        else {'mash': 'not_required', 'minimap2': 'not_required'}
    )
    merge_paf = build_pairwise_overlap_paf(
        fasta_path=args.p_ctg,
        output_directory=str(output_root),
        min_contig_length=args.min_contig_length,
        min_distance=args.min_distance,
        threads_per_alignment=args.threads,
        max_processes=args.process,
        cpu_budget=args.cpu_budget,
        tool_versions=tuple(tool_versions.items()),
    )

    graph_directory = output_root / '04.remove.redundancy'
    graph_directory.mkdir(parents=True, exist_ok=True)
    run_mt2t_overlap_workflow(
        fasta_path=args.p_ctg,
        paf_path=str(merge_paf),
        output_directory=str(graph_directory),
        parameters=parameters,
        tool_versions=tool_versions,
    )


def _tool_version(tool: str) -> str:
    completed = run_command(
        [tool, '--version'],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    output = completed.stdout
    if isinstance(output, bytes):
        return output.decode('utf-8', errors='replace').strip().splitlines()[0]
    return str(output).strip().splitlines()[0]


if __name__ == '__main__':
    main()
