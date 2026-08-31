#!/usr/bin/env bash
# Full PHap v2 phase_reads validation for the current autotetraploid potato data.
# It runs: assignment -> extraction -> hifiasm assembly -> HapHiC scaffolding.
# Run this only after run_cluster_full_potato_v2.sh has completed.
#
# Usage:
#   bash examples/run_phase_reads_full_potato_v2.sh [cluster_output_directory] [output_directory]
#
# Activate the environment containing hifiasm, bwa, samtools, samblaster,
# filter_bam, haphic, seqkit, pigz, and gawk before starting this script.

set -euo pipefail

workspace=/home/pxxiao/test/PHap.cluster
code_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cluster_output=${1:-"${workspace}/02.cluster.v2.full_validation_20260831"}
output_dir=${2:-"${workspace}/03.phase_reads.v2.full_validation_20260831"}
temp_dir="${output_dir}.tmp"

group_file="${cluster_output}/05.rescue/group.reassignment.cluster.txt"
contig_type="${workspace}/contig_depth.txt"
bam_hifi="${workspace}/01.hifi.data/HiFi.p_utg.sort.bam"
bam_ont="${workspace}/01.hifi.data/ONT.p_utg.sort.bam"
bam_hic="${workspace}/01.hifi.data/HiC.p_utg.sort.bam"

# The exact raw inputs used in the previous real-data validation manifests.
raw_hifi=/home/pxxiao/project/10_potato_poly/00_data/02_check_20240123/PB/potato4.hifi.fastq.gz
raw_ont=/home/pxxiao/project/10_potato_poly/02_asm/01_test/01_hifi_ont/00_data/new_ont/100k/hs.ont.100k.fq.gz
raw_hic1=/home/pxxiao/data/04_potato_/WHXWZB-2023100244A/raw_data/MGI/Hi-C/hua-4-4/257.25G/clean/hua-4-4_combine_1.fq.clean.gz
raw_hic2=/home/pxxiao/data/04_potato_/WHXWZB-2023100244A/raw_data/MGI/Hi-C/hua-4-4/257.25G/clean/hua-4-4_combine_2.fq.clean.gz

require_file() {
    [[ -r "$1" ]] || { echo "Missing or unreadable input: $1" >&2; exit 2; }
}

for input in "$group_file" "$contig_type" "$bam_hifi" "$bam_ont" "$bam_hic" \
             "$raw_hifi" "$raw_ont" "$raw_hic1" "$raw_hic2"; do
    require_file "$input"
done

for tool in hifiasm bwa samtools samblaster filter_bam haphic seqkit pigz gawk; do
    command -v "$tool" >/dev/null || {
        echo "Required executable is not on PATH: $tool" >&2
        exit 2
    }
done

# Do not overwrite or silently mix an earlier run. The pipeline itself supports
# --resume, but a new complete validation must begin in an empty directory.
if [[ -e "$output_dir" ]] && [[ -n "$(find "$output_dir" -mindepth 1 -print -quit)" ]]; then
    echo "Output directory is not empty: $output_dir" >&2
    exit 2
fi
mkdir -p "$output_dir" "$temp_dir"

python "${code_root}/PHap.py" phase_reads \
    --bam-hifi "$bam_hifi" \
    --bam-ont "$bam_ont" \
    --bam-hic "$bam_hic" \
    --hifi "$raw_hifi" \
    --ont "$raw_ont" \
    --hic1 "$raw_hic1" \
    --hic2 "$raw_hic2" \
    --contig-type "$contig_type" \
    --group "$group_file" \
    --data-types hifi ont hic \
    --output-dir "$output_dir" \
    --temp-dir "$temp_dir" \
    --stop-after scaffold \
    --seed 100 \
    --progress-every 1000000 \
    --collapsed-policy balanced \
    --unknown-contig-type-policy group \
    --min-group-margin 0.10 \
    --balanced-score-tolerance 0.02 \
    --hifi-min-alignment-length 1000 \
    --hifi-min-identity 0.95 \
    --ont-min-alignment-length 1000 \
    --ont-min-identity 0.75 \
    --hic-min-alignment-length 50 \
    --hic-min-identity 0.90 \
    --hic-assignment-backend fast-v1 \
    --ont-length 1 \
    --ont-quality 0 \
    --extract-backend seqkit \
    --extract-threads 4 \
    --jobs 4 \
    --threads-per-job 10 \
    --haphic-processes 1 \
    --haphic-nx 100 \
    --haphic-nm 3

scaffolds=$(find "${output_dir}/04.scaffold" -type f -path '*/04.build/scaffolds.fa' | wc -l)
checkpoints=$(find "${output_dir}/04.scaffold/.checkpoints" -type f -name '*.complete.json' | wc -l)
[[ "$scaffolds" -eq 48 && "$checkpoints" -eq 48 ]] || {
    echo "Expected 48 completed scaffolds/checkpoints; found scaffolds=${scaffolds}, checkpoints=${checkpoints}" >&2
    exit 1
}

echo "Full phase_reads validation complete: $output_dir"
