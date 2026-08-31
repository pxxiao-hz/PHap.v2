#!/usr/bin/env bash
# Full PHap v2 cluster validation for the current autotetraploid potato data.
# It runs: allelic table -> chromosome extraction -> cluster -> recluster -> rescue.
# Usage: bash examples/run_cluster_full_potato_v2.sh [output_directory]

set -euo pipefail

workspace=/home/pxxiao/test/PHap.cluster
code_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
output_dir=${1:-"${workspace}/02.cluster.v2.full_validation_20260831"}

putg="${workspace}/hs.100k.asm.bp.p_utg.gfa.remove.plastid.contamination.fa"
mt2t="${workspace}/hs.mT2T.v2.fa"
contig_type="${workspace}/contig_depth.txt"
raw_paf="${workspace}/02.cluster/01.putg_vs_mT2T/p_utg_vs_mT2T.paf"
noseq_gfa="${workspace}/hs.100k.asm.bp.p_utg.noseq.gfa"
full_links="${workspace}/full_links.pkl"

require_file() {
    [[ -r "$1" ]] || { echo "Missing or unreadable input: $1" >&2; exit 2; }
}

for input in "$putg" "$mt2t" "$contig_type" "$raw_paf" "$noseq_gfa" "$full_links"; do
    require_file "$input"
done

# Never mix a new test with an old partial run. Choose a new output directory,
# or remove/move a previous test directory yourself after reviewing it.
if [[ -e "$output_dir" ]] && [[ -n "$(find "$output_dir" -mindepth 1 -print -quit)" ]]; then
    echo "Output directory is not empty: $output_dir" >&2
    exit 2
fi
mkdir -p "$output_dir"

python "${code_root}/PHap.py" cluster \
    --p_utg "$putg" \
    --mT2T "$mt2t" \
    --contig_type "$contig_type" \
    --paf "$raw_paf" \
    --gfa "$noseq_gfa" \
    --full_links "$full_links" \
    --output_dir "$output_dir" \
    --stop_after rescue \
    --threads 4 \
    --chr_num 12 \
    --top_n 4 \
    --hic_link_normalization dosage \
    --cluster_balance_weight 1.0 \
    --cluster_refinement_rounds 10 \
    --cluster_phase_block_rounds 10 \
    --cluster_min_phase_block_gain 0.005 \
    --cluster_max_phase_boundary_relaxations 1 \
    --cluster_max_phase_boundary_overlap 0.30 \
    --cluster_max_constraint_relaxation_rounds 50 \
    --cluster_min_constraint_relaxation_gain 0.0 \
    --cluster_min_constraint_relaxation_links 5.0 \
    --cluster_min_constraint_relaxation_margin 0.10 \
    --cluster_max_constraint_relaxations_per_move 2 \
    --cluster_max_constraint_relaxation_overlap 0.30 \
    --cluster_max_phase_interval_rounds 4 \
    --cluster_max_phase_interval_relaxations 2 \
    --cluster_max_backtracks 1000000 \
    --cluster_constraint_relaxation weakest \
    --cluster_protected_long_unitig_length 5000000 \
    --cluster_protected_length_ratio 5.0 \
    --cluster_protected_short_overlap 0.50 \
    --cluster_protected_direct_overlap_bp 1000000 \
    --cluster_protected_direct_short_overlap 0.20 \
    --recluster_min_adjusted_links 5.0 \
    --recluster_min_group_margin 0.10 \
    --recluster_min_allelic_block_anchors 1 \
    --recluster_max_allelic_block_configurations 256 \
    --recluster_max_allelic_block_movable_length 5000000 \
    --recluster_min_group_bp_ratio 0.25 \
    --recluster_min_assigned_fraction 0.0 \
    --recluster_max_rounds 10 \
    --recluster_refinement_rounds 4 \
    --recluster_low_confidence_policy defer \
    --recluster_unknown_dosage_policy haplotig \
    --recluster_seed_review weak \
    --recluster_trusted_seed_bases hic_supported \
    --recluster_reviewed_seed_fallback retain \
    --rescue_min_adjusted_links 5.0 \
    --rescue_min_chromosome_margin 0.10 \
    --rescue_min_group_margin 0.10 \
    --rescue_min_assigned_fraction 0.0 \
    --rescue_max_rounds 10 \
    --rescue_low_confidence_policy defer \
    --rescue_unknown_dosage_policy defer

group_file="${output_dir}/05.rescue/group.reassignment.cluster.txt"
[[ -s "$group_file" ]] || { echo "Missing final group file: $group_file" >&2; exit 1; }

group_fastas=$(find "${output_dir}/05.rescue" -maxdepth 1 -type f \
    -name 'chr*_group*.reassignment.fa' | wc -l)
[[ "$group_fastas" -eq 48 ]] || {
    echo "Expected 48 final group FASTAs, found ${group_fastas}" >&2
    exit 1
}

echo "Cluster validation complete: $output_dir"
echo "Final phase_reads input: $group_file"
