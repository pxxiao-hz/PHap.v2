# PHap v2 Constraint-Aware Clustering

## Purpose

`03.cluster` assigns every unitig in the accepted allelic table to exactly the
number of haplotype groups implied by its dosage. Allelic-table pairs are
mutual-exclusion constraints: two enforced allelic unitigs may never share a
group. A unitig is not discarded because it has weak or absent Hi-C links.

## Algorithm

1. Read each chromosome FASTA once and retain only table unitigs for the seed
   clustering stage.
2. Convert contig types to exact group dosage: haplotig 1, diplotig 2,
   triplotig 3, and tetraplotig 4.
3. Build a conflict graph from unitigs that co-occur in an allelic-table row.
4. Detect conflict cliques whose summed dosage exceeds ploidy. In the default
   `weakest` mode, remove only the least-supported edge in each impossible
   clique, ranked first by the smaller directional table-overlap ratio and then
   by overlap bp. Substantial direct-projection edges (at least 1 Mb and at
   least 20% of the shorter unitig by default) are protected and cannot be
   removed. Envelope-only edges remain eligible for relaxation. Record every
   removed edge. `fail` mode stops instead.
5. By default, divide each Hi-C count by the product of the two unitig dosages,
   preventing collapsed unitigs from receiving an artificial interaction
   advantage. `--hic_link_normalization raw` keeps the original read-pair
   counts for a controlled comparison.
6. Select a full-dosage table row with strong external Hi-C as an initial
   anchor, then solve the exact-dosage graph-coloring problem with forward
   checking and deterministic DSATUR-style variable ordering.
7. Refine legal assignments using normalized Hi-C density and group bp load.
8. Apply Kempe-component swaps. These exchange two group labels across an
   entire swappable conflict component, preserving every enforced conflict
   while escaping single-unit local optima.
9. Review assignments that are strongly contradicted by normalized Hi-C. A
   candidate must pass minimum adjusted-link and group-margin thresholds,
   increase Hi-C cohesion, and not decrease the global objective. It may release only a
   bounded number of non-protected table edges whose overlap fraction is below
   the configured ceiling. Every released edge and group change is audited.
10. Scan allelic-table coordinates for phase-block boundaries. At each boundary,
   test every pairwise suffix group permutation and accept only moves that
   improve the internal Hi-C/balance objective by the configured minimum. A
   move may relax at most the configured number of non-protected boundary
   edges, and only when their table overlap fraction is below the configured
   ceiling. Unitigs always move as complete records; they are never split at a
   boundary. Record every move and relaxed edge.
11. Test pairwise group permutations inside intervals with sequence on both
   sides. This two-boundary search can repair an internal phase block without
   changing an already-correct distal suffix. Coordinate buckets unaffected by
   a group-pair swap are skipped without changing its score.
12. Validate dosage, enforced conflicts, and assignment completeness before
   replacing output files.

The global refinement objective is:

```text
Hi-C cohesion - balance_weight * mean squared relative group-bp deviation
```

Both terms are dimensionless. The default `balance_weight=1.0` worked better
than 0.10 on the current data without forcing the four groups to equal length
when one haplotype is represented by an exceptionally long unitig.

## Main Parameters

| Public option | Default | Meaning |
| --- | ---: | --- |
| `--hic_link_normalization` | `dosage` | Hi-C counts used by stages 03-05: dosage-corrected or `raw` |
| `--cluster_balance_weight` | 1.0 | Trade-off between Hi-C cohesion and seed bp balance |
| `--cluster_refinement_rounds` | 10 | Maximum rounds for each refinement stage |
| `--cluster_phase_block_rounds` | 10 | Maximum accepted-boundary search rounds |
| `--cluster_min_phase_block_gain` | 0.005 | Minimum objective gain for a block move |
| `--cluster_max_phase_boundary_relaxations` | 1 | Maximum weak boundary edges relaxed per move |
| `--cluster_max_phase_boundary_overlap` | 0.30 | Maximum endpoint overlap fraction for a boundary relaxation |
| `--cluster_max_constraint_relaxation_rounds` | 50 | Maximum evidence-directed weak-edge rounds; stops early at stability |
| `--cluster_min_constraint_relaxation_gain` | 0.0 | Minimum global objective gain for a weak-edge move |
| `--cluster_min_constraint_relaxation_links` | 5.0 | Minimum dosage-adjusted links supporting candidate groups |
| `--cluster_min_constraint_relaxation_margin` | 0.10 | Minimum normalized density margin for candidate groups |
| `--cluster_max_constraint_relaxations_per_move` | 2 | Maximum weak edges released by one unitig move |
| `--cluster_max_constraint_relaxation_overlap` | 0.30 | Maximum endpoint overlap fraction for a released edge |
| `--cluster_max_phase_interval_rounds` | 4 | Maximum internal-interval search rounds |
| `--cluster_max_phase_interval_relaxations` | 2 | Maximum weak edges released by one interval move |
| `--cluster_max_backtracks` | 1,000,000 | Hard limit for exact constraint search |
| `--cluster_constraint_relaxation` | `weakest` | Resolve impossible cliques or stop with `fail` |
| `--cluster_protected_direct_overlap_bp` | 1,000,000 | Minimum accepted-block overlap for a protected direct edge |
| `--cluster_protected_direct_short_overlap` | 0.20 | Minimum direct overlap fraction of the shorter unitig |

The clustering utility also exposes phase-block controls directly:

| Utility option | Default | Meaning |
| --- | ---: | --- |
| `--hic-link-normalization` | `dosage` | Divide by endpoint dosages or retain raw counts |
| `--max-phase-block-rounds` | 10 | Maximum accepted-boundary search rounds |
| `--min-phase-block-gain` | 0.005 | Minimum dimensionless objective improvement |
| `--max-phase-boundary-relaxations` | 1 | Maximum weak edges relaxed by one block move |
| `--max-phase-boundary-overlap` | 0.30 | Maximum overlap fraction on either endpoint of a relaxable edge |
| `--max-constraint-relaxation-rounds` | 50 | Maximum evidence-directed weak-edge rounds |
| `--min-constraint-relaxation-gain` | 0.0 | Minimum global objective gain |
| `--min-constraint-relaxation-links` | 5.0 | Minimum candidate adjusted links |
| `--min-constraint-relaxation-margin` | 0.10 | Minimum candidate density margin |
| `--max-constraint-relaxations-per-move` | 2 | Maximum weak edges released by one move |
| `--max-constraint-relaxation-overlap` | 0.30 | Maximum endpoint overlap fraction for those edges |
| `--max-phase-interval-rounds` | 4 | Maximum internal-interval rounds |
| `--max-phase-interval-relaxations` | 2 | Maximum weak edges released by one interval move |
| `--allelic-pairs` | none | Pair-evidence sidecar generated with the allelic table |
| `--protected-direct-overlap-bp` | 1,000,000 | Minimum direct-projection overlap for protection |
| `--protected-direct-short-overlap` | 0.20 | Minimum direct overlap fraction of the shorter unitig |

## Outputs Per Chromosome

| File | Meaning |
| --- | --- |
| `g1.fa` ... `g4.fa` | Legacy-compatible group FASTA files |
| `g1.txt` ... `g4.txt` | Legacy-compatible unitig lists |
| `group.cluster.txt` | Legacy group summary |
| `g1g2g3g4.txt` / `.fa` | Legacy combined outputs |
| `cluster_assignments.tsv` | Per-unitig dosage, groups, Hi-C support, margin, and constraint degree |
| `cluster_summary.json` | Inputs, parameters, group sizes, search/refinement metrics, and validation |
| `cluster_constraint_violations.tsv` | Must contain only its header in an accepted run |
| `cluster_relaxed_constraints.tsv` | Every deliberately relaxed edge and its evidence |
| `cluster_protected_constraints.tsv` | Every non-relaxable edge, protection reason, table overlap, and direct-projection support |
| `cluster_phase_block_moves.tsv` | Boundary, group permutation, objective change, and relaxed edges for every accepted phase-block move |
| `cluster_constraint_relaxation_moves.tsv` | Unitig move, local Hi-C evidence, objective change, and released weak edges |
| `cluster_phase_interval_moves.tsv` | Two interval boundaries, group permutation, objective change, and released edges |

## Acceptance Checks

An accepted chromosome has all three summary values equal to zero:

```text
validation.allelic_conflicts
validation.dosage_errors
validation.unassigned_table_unitigs
```

Also review `cluster_relaxed_constraints.tsv`; a nonempty file is an input-table
inconsistency requiring biological review, even though the remaining enforced
problem is valid. Group bp need not be equal when long haplotype-specific
unitigs differ in size or the conservative allelic table covers haplotypes
unevenly. The following `04.recluster` stage adds non-table chromosome unitigs.
Its confidence and completeness rules are documented in
[RECLUSTER_V2.md](RECLUSTER_V2.md).
