#!/bin/bash

set -euo pipefail

# Usage message
usage() {
    echo "Usage: $0 <genome> <hic1> <hic2> [--threads <num_threads>]"
    echo "Arguments:"
    echo "  <genome>: Path to the reference genome"
    echo "  <hic1>: Path to the first Hi-C read file"
    echo "  <hic2>: Path to the second Hi-C read file"
    echo "Options:"
    echo "  --threads <num_threads>: Number of threads for parallel processing (default: 32)"
}

# Default value for threads
num_threads=32

# Check for correct number of arguments
if [ "$#" -lt 3 ]; then
    echo "Error: Incorrect number of arguments!"
    usage
    exit 1
fi

# Assign positional arguments first.
genome="$1"
hic1="$2"
hic2="$3"
shift 3

# Parse options.
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --threads)
            if [ "$#" -lt 2 ]; then
                echo "Error: --threads requires a value." >&2
                exit 1
            fi
            num_threads="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

filter_bam_command="${HAPHIC_FILTER_BAM:-filter_bam}"
for command in bwa samblaster samtools "$filter_bam_command"; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "Error: required executable not found on PATH: $command" >&2
        exit 127
    fi
done

# Index the reference genome
bwa index "$genome"

# Perform read alignment and processing
bwa mem -5SP -t "$num_threads" "$genome" "$hic1" "$hic2" | samblaster | samtools view - -@ "$num_threads" -S -h -b -F 3340 -o HiC.bam
"$filter_bam_command" HiC.bam 1 --nm 3 --threads "$num_threads" | samtools view - -b -@ "$num_threads" -o HiC.filtered.bam
