#!/bin/bash

# Usage message
usage() {
    echo "Usage: $0 <genome> <hic1> <hic2> [--threads <num_threads>]"
    echo "Arguments:"
    echo "  <genome>: Path to the reference genome"
    echo "  <hic1>: Path to the first Hi-C read file"
    echo "  <hic2>: Path to the second Hi-C read file"
    echo "Options:"
    echo "  --threads <num_threads>: Number of threads for parallel processing (default: 32)"
    exit 1
}

# Default value for threads
num_threads=32

# Check for correct number of arguments
if [ "$#" -lt 3 ]; then
    echo "Error: Incorrect number of arguments!"
    usage
fi

# Parse command line options
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --threads)
            num_threads="$2"
            shift
            shift
            ;;
        *)
            break
            ;;
    esac
done

# Assign non-option arguments to variables
genome="$1"
hic1="$2"
hic2="$3"

# Required executables must be available on PATH. FILTER_BAM may be used to
# provide an explicit path to HapHiC's filter_bam utility.
filter_bam_command="${FILTER_BAM:-filter_bam}"
for command in bwa samblaster samtools "$filter_bam_command"; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "Error: required executable not found: $command" >&2
        exit 1
    fi
done

# Index the reference genome
bwa index "$genome"

# Perform read alignment and processing
bwa mem -5SP -t "$num_threads" "$genome" "$hic1" "$hic2" | samblaster | samtools view - -@ "$num_threads" -S -h -b -F 3340 -o HiC.bam
"$filter_bam_command" HiC.bam 1 --nm 3 --threads "$num_threads" | samtools view - -b -@ "$num_threads" -o HiC.filtered.bam
