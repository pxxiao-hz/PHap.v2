# PHap v2 Chromosome-Unassigned Unitig Rescue

## Purpose

`05.rescue` is the final conservative assignment step. It considers only the
unitigs in `02.chr_seq/un_chr.fa`. Assignments and explicit deferrals from
`04.recluster` are immutable.

This scope is important: a unitig deferred by chromosome-level reclustering
must not be reconsidered against all chromosomes, because cross-chromosome
Hi-C noise could otherwise move it away from its sequence-supported
chromosome.

## Algorithm

1. Load all `04.recluster/chr*/recluster_assignments.tsv` files and verify
   their group tables, dosage, FASTA membership, and chromosome ownership.
2. Parse only `un_chr.fa` as rescue candidates. Verify that the union of these
   candidates and all chromosome inputs exactly equals the original p_utg
   FASTA, including sequence lengths.
3. Keep only Hi-C pairs involving a rescue candidate. The default divides
   every count by the product of endpoint dosages; the public
   `--hic_link_normalization raw` option retains raw read-pair counts.
4. For each chromosome, choose exactly `dosage` groups by normalized Hi-C
   density. Compare the best chromosome with the second-best chromosome, then
   compare the weakest selected group with the strongest unselected group
   inside the selected chromosome.
5. Require minimum adjusted links, positive support for every selected group,
   a chromosome margin, and the group-density margin. Accept candidates in
   batches and permit later rounds to use only earlier high-confidence
   assignments.
6. Defer `other`, `replotig`, missing, or unsupported dosage types by default.
   They can be rejected as input errors or explicitly treated as haplotigs.
7. Validate exact dosage, immutable recluster results, one-chromosome group
   membership, and the complete assigned/unassigned partition before replacing
   output files.

## Defaults

| Public option | Default | Meaning |
| --- | ---: | --- |
| `--hic_link_normalization` | `dosage` | Hi-C counts used by stages 03-05: dosage-corrected or `raw` |
| `--rescue_min_adjusted_links` | 5.0 | Minimum adjusted support in selected groups |
| `--rescue_min_chromosome_margin` | 0.10 | Minimum best-vs-second chromosome density margin |
| `--rescue_min_group_margin` | 0.10 | Minimum dosage-boundary group density margin |
| `--rescue_min_assigned_fraction` | 0.0 | Optional selected-link fraction filter |
| `--rescue_max_rounds` | 10 | Maximum high-confidence propagation rounds |
| `--rescue_low_confidence_policy` | `defer` | Retain ambiguity instead of forcing the best group |
| `--rescue_unknown_dosage_policy` | `defer` | Defer unsupported types; alternatives are `error` and `haplotig` |

## Outputs

| File | Meaning |
| --- | --- |
| `group.reassignment.cluster.txt` | Final downstream membership contract for every chromosome group |
| `chr*_group*.reassignment.fa` | Final group FASTA files, copied from trusted step 04 outputs and appended with rescued candidates |
| `rescue_assignments.tsv` | One row per original p_utg unitig, including source, status, dosage, groups, reason, and decision evidence |
| `rescue_summary.json` | Input partition, candidate decisions, group sizes, Hi-C counts, and validation totals |
| `unassigned_unitigs.txt` / `.fa` | Step 04 deferrals plus rescue candidates that remain unassigned |
| `rescue_validation.tsv` | Must contain only its header in an accepted run |
| `split_clms/chr*_group*.clm` | Final CLM files generated from one scan of the source CLM |
| `clm_split.summary.json` | Input/output CLM record counts for all final groups |

The obsolete combined multi-group FASTA is not generated. It duplicated every
collapsed sequence and substantially increased both peak memory and disk use.

The standalone utility uses `--hic-link-normalization`. Legacy TSV column
names containing `adjusted` are retained for compatibility; their values use
the selected normalization, which is recorded in `rescue_summary.json`.

## Standalone Command

```shell
python PHap.v2/utils/unchr_recluster.py \
  --assembly-fasta p_utg.fa \
  --candidate-fasta 02.cluster.v2/02.chr_seq/un_chr.fa \
  --contig-type contig_depth.txt \
  --full-links full_links.pkl \
  --recluster-dir 02.cluster.v2/04.recluster \
  --output-dir 02.cluster.v2/05.rescue \
  --ploidy 4
```

The public `phap cluster` workflow runs the same command and then invokes the
single-pass CLM splitter automatically.

## Acceptance Checks

- `validation.violations`, dosage, changed-step-04, partition, and length
  errors must all be zero.
- Every assigned row must have `group_count == dosage`, and all listed groups
  must belong to one chromosome.
- Every original p_utg ID must occur either in final groups or in
  `unassigned_unitigs.fa`, never both.
- Review low-support, group-margin, unsupported-type, and no-Hi-C counts before
  changing the default `defer` policies.
