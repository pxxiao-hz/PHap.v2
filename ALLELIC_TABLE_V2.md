# PHap v2 Allelic Unitig Table

## Goal

The allelic table is treated as a set of high-confidence mutual-exclusion
constraints, not as a way to force every unitig into a row. An uncertain
projection is rejected and reported. This deliberately trades some recall for
precision because a false allelic edge can corrupt all downstream haplotype
groups.

## Algorithm

1. Filter raw PAF records by primary status, alignment length, identity, and
   MAPQ. Secondary repeat alternatives do not vote for chromosome assignment.
2. Aggregate records by query and reference chromosome, not by strand. The
   chromosome score integrates the strongest identity/MAPQ-weighted evidence
   at each query position, so overlapping records and opposite strands are not
   double counted. No global monotonicity is assumed.
3. Select one best chromosome per unitig using cumulative local evidence and
   the best-versus-second chromosome margin. A unitig may contain both strands,
   increasing and decreasing runs, inversions, and separated mappings. Query
   coverage and target-span coverage remain QC metrics but do not reject it.
4. Write every reliable local record on the selected chromosome as
   `segmented` evidence. Nearby target intervals are grouped locally. A group
   below the projection-density threshold is recursively split at its largest
   unsupported reference gaps rather than discarded or bridged. Mixed strand
   and nonmonotonic order are valid for `segmented` blocks. Identity, MAPQ,
   minimum aligned length, and actual reference support are still enforced.
   Blocks above the adaptive anchor-length threshold are labeled `anchor`;
   smaller accepted blocks are `supporting`.
5. Sweep all projection endpoints to create atomic half-open reference
   intervals. Two unitigs occur in the same row only where their projected
   intervals actually overlap. Adjacent intervals never form an allelic pair.
6. When a hifiasm GFA is supplied, reject an atomic interval if any proposed
   allelic pair has a direct `L` edge. A direct assembly-graph edge represents
   sequence adjacency/overlap, so it conflicts with treating the pair as
   mutually exclusive alleles at that interval.
7. Sum dosage in every interval. When dosage exceeds ploidy, retain a subset
   only if the best dosage-valid solution exceeds the second-best solution by
   the configured confidence margin. Otherwise omit the interval.
8. Merge adjacent intervals only when their selected unitig set is identical.

Hi-C is intentionally not used to create the table. It remains independent
evidence for clustering and validation, avoiding circular inference.

## Conservative Defaults

| Check | Default | Purpose |
| --- | ---: | --- |
| Alignment length | 1 kb | Remove tiny local hits |
| Alignment identity | 0.90 | Minimum identity for each local record |
| Alignment MAPQ | 20 | Remove ambiguous repeat placements |
| Query length | 20 kb | Avoid weak short-unitig placement |
| Target evidence length | 20 kb | Require enough local evidence to assign a chromosome |
| Aggregate identity | 0.90 | Do not reimpose a stricter whole-unitig divergence filter |
| Best-reference margin | 0.05 | Reject multi-reference ambiguity |
| Projection gap | 20 kb | Initial grouping distance for local target blocks |
| Projection blocks | 10 | Flag highly fragmented unitigs for QC |
| Projection aligned length | 10 kb | Remove tiny projection blocks |
| Projection coverage | 0.70 | Split sparse groups until local projections are dense |
| Block identity | 0.90 | Reject low-identity blocks independently |
| Block mean MAPQ | 20 | Reject weakly placed blocks independently |
| Block collinearity | 0.80 | QC/legacy-mode filter; not a segmented-placement requirement |
| Anchor block length | max(100 kb, 0.5% query) | Label major projection evidence |
| Atomic interval length | 10 kb | Remove tiny shared intervals |
| Dosage resolution margin | 0.10 | Avoid guessing over-capacity subsets |

For a stricter precision run, start with:

```shell
python PHap.py allelic_table \
  --p_utg p_utg.fa --mT2T mT2T.fa --contig_type contig_depth.txt \
  --paf p_utg_vs_mT2T.paf \
  --min-alignment-mapq 30 \
  --min-block-identity 0.95 \
  --min-segment-length 20000 \
  --min-resolution-margin 0.20
```

Thresholds should be calibrated against known haplotypes, parental markers,
or simulated truth. No parameter set can establish biological allelism from
sequence projection alone.

## Outputs

| File | Meaning |
| --- | --- |
| `corrected_allelic_table.txt` | Headerless table consumed by clustering |
| `allelic.ctg.table.v2` | Compatibility alias of the primary table |
| `putg_vs_mT2T.collinear.paf` | Compatibility filename containing all reliable segmented records on one selected chromosome; records are not required to be globally collinear |
| `collinear_chain.qc.tsv` | Chromosome-selection decision, both-strand metrics, switches, margin, and rejection reason |
| `chromosome_selection.summary.json` | Active selection parameters, accepted sequence/record counts, mixed-strand counts, and rejection reasons |
| `unitig_projections.tsv` | Localized block-level query/target coverage, contribution, identity, MAPQ, strand, collinearity, class, and confidence |
| `allelic_table.qc.tsv` | Every atomic interval, dosage, decision, and margin |
| `rejected_projections.tsv` | Projection rejection reason and detail |
| `allelic_pairs.tsv` | Pair overlap bp and directional overlap ratios |
| `allelic_pairs.gfa_validation.tsv` | Per-pair graph topology and read-assignment evidence |
| `allelic_pairs.gfa_validation.summary.json` | Counts for each GFA evidence class |
| `allelic_table.summary.json` | Counts and table-generation parameters |
| `run_manifest.json` | Input fingerprints, parameters, and output paths |

## Acceptance Checks Before Clustering

1. `dosage_gt_ploidy` must be zero in the final table.
2. Every row must have positive length and unique unitig IDs.
3. Every pair in `allelic_pairs.tsv` must have positive `overlap_bp`.
4. Inspect rejection counts by reason; a large unknown-dosage count means the
   dosage input needs correction, not that the unitigs should be guessed as
   haplotigs.
5. Review low directional overlap ratios and high-confidence Hi-C links. These
   are candidates for structural variation, repeats, chimerism, or an incorrect
   mT2T projection.
6. Compare table precision with parental markers or known phased sequence
   before enabling the complete clustering workflow.

## Interpreting GFA Evidence

The GFA check is asymmetric: some graph structures can reject a proposed
allelic relation with high confidence, while absence of such a structure does
not prove the relation.

| Classification | Interpretation |
| --- | --- |
| `conflict_direct_link` | The two unitigs have a direct `L` edge and should not be treated as allelic over that graph overlap |
| `support_bubble_like` | The unitigs share at least two graph neighbors; useful positive bubble-like evidence, but still review repeats and complex graph regions |
| `review_shared_reads` | The unitigs share hifiasm `A`-record reads without stronger bubble evidence; inspect before use as proof |
| `graph_neutral` | The GFA neither supports nor contradicts the pair; rely on collinearity, dosage, markers, and independent Hi-C evidence |

A no-sequence GFA is sufficient because this validation uses `S`, `L`, and
`A` records rather than segment sequence. GFA validation must not be described
as a biological ground truth set.

## Clustering Changes in v2

The v2 cluster entry point consumes the new table directly and can stop after
any completed stage with `--stop_after`; use `--stop_after cluster` to avoid
running recluster and rescue. `--table_only` remains a compatibility alias for
stopping after table generation. The rewritten constraint solver assigns
the exact dosage, retains table unitigs without Hi-C, and never overrides an
enforced table constraint. Globally impossible conflict cliques either stop the
run or cause one weakest edge to be relaxed and explicitly reported. See
[CLUSTER_V2.md](CLUSTER_V2.md) for the algorithm and QC outputs.
