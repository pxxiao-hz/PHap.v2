# PHap v2 Chromosome Reassignment

## Purpose

`04.recluster` reviews uncertain `03.cluster` assignments and adds
chromosome-assigned unitigs that were not used as high-confidence seeds.
Only seeds whose cluster evidence passes the configured Hi-C thresholds are
immutable. Other seeds are removed from the propagation anchors and rescored
without allelic-table constraints.

Uncertain unitigs are retained in explicit unassigned FASTA and TSV outputs.
They are not silently deleted and are not guessed into a haplotype by default.

## Algorithm

1. Validate every seed ID, its dosage, and its matching row in
   `03.cluster/cluster_assignments.tsv`.
2. Keep only trusted `hic_supported` seeds whose adjusted links, assigned-link
   fraction, and group margin pass the recluster thresholds as fixed anchors.
3. Review all other seeds using chromosome Hi-C alone. A strong contradictory
   result changes the groups; an inconclusive result retains the original
   groups by default but is not allowed to propagate uncertainty.
4. Convert haplotig, diplotig, triplotig, and tetraplotig types to dosages one
   through four. Missing types can stop the run or be explicitly recorded as
   inferred haplotigs.
5. Load only Hi-C pairs whose endpoints occur in the chromosome FASTA. The
   default divides each count by the product of endpoint dosages; the public
   `--hic_link_normalization raw` option retains raw read-pair counts.
6. Sum adjusted links to each current group and normalize by that group's
   restriction-site count.
7. Select exactly `dosage` groups. Require minimum adjusted links, positive
   support for every selected group, and the configured density margin between
   the weakest selected and strongest unselected group.
8. Accept candidates in batches, then repeat. A later round may use only the
   assignments accepted by earlier rounds, avoiding length-order dependence.
9. Synchronously rescore every non-fixed unitig after propagation. Repeat four
   times (equivalent to complete recluster runs 2 through 5), stopping early
   when no group changes. Repeated states are detected and reported as an
   oscillation instead of being silently accepted.
10. Assign tetraplotigs to all four groups by dosage even without Hi-C.
11. Validate fixed-seed preservation, exact dosage, unique memberships, and the
   assigned/unassigned input partition before replacing outputs.
12. Read chromosome-local allelic-table rows after Hi-C refinement. Treat
   trusted assignments as anchors and jointly enumerate dosage-exact,
   mutually-exclusive group configurations for weak members. Accept only a
   unique constraint-forced completion or a block-margin-supported best
   configuration. Pairs relaxed by `03.cluster` are excluded so a weak
   table edge is not silently restored.

## Defaults

| Public option | Default | Meaning |
| --- | ---: | --- |
| `--hic_link_normalization` | `dosage` | Hi-C counts used by stages 03-05: dosage-corrected or `raw` |
| `--recluster_min_adjusted_links` | 5.0 | Minimum adjusted support in selected groups |
| `--recluster_min_group_margin` | 0.10 | Minimum density boundary margin |
| `--recluster_min_allelic_block_anchors` | 1 | Minimum anchored members for joint block resolution |
| `--recluster_max_allelic_block_configurations` | 256 | Safety cap on exact block enumeration |
| `--recluster_max_allelic_block_movable_length` | 5000000 | Prevent a chromosome-local allelic block from moving an entire unitig at least this long; `0` disables the protection |
| `--recluster_min_group_bp_ratio` | 0.25 | Reject a chromosome if its smallest haplotype group is below this fraction of median group bp; `0` disables the safety check |
| `--recluster_min_assigned_fraction` | 0.0 | Optional absolute selected-link fraction filter |
| `--recluster_max_rounds` | 10 | Maximum high-confidence propagation rounds |
| `--recluster_refinement_rounds` | 4 | Stability rounds after initial propagation; rounds 2 through 5 |
| `--recluster_low_confidence_policy` | `defer` | Defer ambiguity or force best available groups |
| `--recluster_unknown_dosage_policy` | `haplotig` | Auditably infer missing types or stop with `error` |
| `--recluster_seed_review` | `weak` | Review weak/constraint-driven cluster seeds; `off` restores immutable legacy seeds |
| `--recluster_trusted_seed_bases` | `hic_supported` | Cluster evidence bases eligible to become fixed anchors |
| `--recluster_reviewed_seed_fallback` | `retain` | Retain an original seed assignment when Hi-C review is inconclusive; `defer` removes it from final groups |

The selected-link fraction is not filtered by default. A valid haplotig can
have background links to all four groups while retaining a clear normalized
density advantage for one group.

The standalone utility uses `--hic-link-normalization`. Legacy TSV column
names containing `adjusted` are retained for compatibility; their values use
the selected normalization, which is recorded in `recluster_summary.json`.

## Outputs Per Chromosome

| File | Meaning |
| --- | --- |
| `group1.reassignment.fa` ... | Legacy-compatible assigned group FASTA |
| `group.reassignment.cluster.txt` | Downstream group membership contract |
| `recluster_assignments.tsv` | Every input unitig, original seed groups, cluster basis, review outcome, dosage, decision evidence, and final evidence |
| `recluster_refinement.tsv` | Every group change made in stability rounds, with old/new groups and Hi-C evidence |
| `allelic_block_reassignments.tsv` | Every accepted joint block, anchors, alternatives, acceptance basis, and changed groups |
| `allelic_block_protections.tsv` | Blocks where chromosome-scale unitigs were retained as anchors instead of being moved by local block evidence |
| `recluster_summary.json` | Input/assigned/unassigned bp, rounds, reasons, groups, and validation |
| `unassigned_unitigs.txt` / `.fa` | Deferred sequence retained with no guessed haplotype |
| `recluster_validation.tsv` | Must contain only its header in an accepted run |

## Acceptance Checks

- `validation.violations`, `dosage_errors`, `changed_fixed_seeds`,
  `group_bp_imbalances`, and
  `partition_errors` must all be zero.
- Every assigned row must have `group_count == dosage`.
- Review all `low_group_margin` unitigs before using
  `--recluster_low_confidence_policy best`.
- Review `reviewed_seed.changed_unitigs`; these are explicit corrections to
  the `03.cluster` result, not newly added unitigs.
- `refinement.stable` must be true and `refinement.oscillation` must be false.
- Missing dosage records must be reviewed through `dosage_source` and
  `inferred_haplotig_unitigs`.
- Group bp is a QC signal, not an assignment rule. Strong Hi-C is not
  overridden merely to equalize FASTA sizes.
