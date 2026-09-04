# PHap v2 general ploidy support

## Scope

The v2 workflow accepts any integer ploidy of at least two.  `--top_n` is the
expected haplotype count, while `--chr_num` is the chromosome count in one
basic genome set.  A triploid species with 12 basic chromosomes therefore uses
`--top_n 3 --chr_num 12` and produces 36 final chromosome groups.

The allelic-table, clustering, chromosome reassignment, rescue, and read
phasing stages all use the configured or inferred group count. Synthetic tests
cover triploid and pentaploid clustering, dosage-five allelic-table input, and
the complete dosage-six path through these stages. Biological end-to-end
calibration is still limited to autotetraploid potato; runs at another ploidy
or in another species remain
experimental until their QC and independent validation are reviewed.

## Dosage input

The `contig_type` column accepts the historical names:

| Dosage | Label |
| ---: | --- |
| 1 | `haplotig` |
| 2 | `diplotig` |
| 3 | `triplotig` |
| 4 | `tetraplotig` |
| 5 | `pentaplotig` |
| 6 | `hexaplotig` |

For any positive dosage, a bare integer or generic labels `dosage_N`,
`copy_N`, and `Nx` are also accepted, for example `7`, `dosage_7`, `copy_8`,
or `9x`.  A parsed dosage
greater than the configured ploidy is rejected rather than silently converted
to a haplotig.

Classify contigs from depth using a sample-specific single-copy peak:

```shell
python shell/dosage.analysis.contig.type.identified.py \
  --input_file aln.sort.clean.pandepth.win.stat.txt \
  --pandepth \
  --ploidy 6 \
  --base-depth 21
```

The classifier assigns the nearest integer dosage using half-copy boundaries.
Depth below 0.5 times the single-copy peak is `other`; depth above the
configured ploidy plus half a copy is `replotig`.  Inspect the empirical depth
histogram before choosing `--base-depth`; an incorrect peak directly corrupts
all downstream dosage constraints.

## No-collapse mode

Use `--no-collapse` when the input assembly is known to contain no collapsed
unitigs. `--contig_type` is then omitted, and every unitig is assigned dosage
one throughout allelic-table generation, clustering, chromosome reassignment,
rescue, and read phasing. Consequently, an assigned unitig can occur in
exactly one haplotype group. `--no-collapse` and `--contig_type` are mutually
exclusive.

```shell
python PHap.py cluster \
  --p_utg p_utg.fa \
  --mT2T mT2T.fa \
  --full_links full_links.pkl \
  --top_n 3 \
  --chr_num 12 \
  --no-collapse \
  --output_dir 02.cluster.triploid.no_collapse
```

This is Hi-C-driven single-copy clustering with the existing mT2T chromosome
placement and allelic mutual-exclusion constraints. Because every dosage is
one, dosage-normalized Hi-C counts equal raw counts.

The public workflow still runs `04.recluster` and `05.rescue`. In this mode
they do not perform multi-group collapsed-unitig recovery: step 04 assigns
chromosome-local unitigs omitted from the seed table, and step 05 attempts to
place unitigs that lacked sequence-based chromosome placement. Keeping them
prevents those sequences from disappearing from the final partition. Use
`--stop_after cluster` only when an intentionally incomplete, seed-table-only
result is desired.

Read assignment must use the same contract:

```shell
python PHap.py phase_reads \
  --group 02.cluster.triploid.no_collapse/05.rescue/group.reassignment.cluster.txt \
  --no-collapse \
  ...
```

The group file is rejected if any unitig occurs in more than one group. The
collapsed-read balancing policy is therefore unused. If the assembly actually
contains collapsed or duplicated sequence, forcing this mode can assign it to
only one haplotype and bias both group content and read depth.

## Workflow

```shell
python PHap.py cluster \
  --p_utg p_utg.fa \
  --mT2T mT2T.fa \
  --contig_type contig_depth.txt \
  --full_links full_links.pkl \
  --top_n 6 \
  --chr_num 12 \
  --output_dir 02.cluster.hexaploid
```

For non-tetraploids the combined step-03 outputs are `all_groups.txt` and
`all_groups.fa`; per-group files remain `g1` through `gP`.  Step 04 similarly
uses `all_groups.reassignment.fa`.  The final group contract always uses
`chrN_group1` through `chrN_groupP`, which lets `phase_reads` infer ploidy
without a separate option.

## Required validation

Start with one chromosome and inspect:

- rejected projections and dosage-over-ploidy records;
- exact `dosage == group_count` in cluster, recluster, and rescue tables;
- constraint violations and relaxed allelic edges;
- per-group bp balance and Hi-C density margins;
- deferred or unsupported contigs;
- read assignment conservation, exclusivity, and estimated group depth.

Higher ploidy enlarges a dosage-`d` unitig domain to `C(ploidy, d)`.  Increase
`--cluster_max_backtracks` or `--recluster_max_allelic_block_configurations`
only after inspecting the component that reached the limit; raising a limit
does not resolve contradictory allelic evidence.  Hi-C and balance thresholds
were calibrated on potato and may require species-specific sensitivity tests.

PHap is designed for autopolyploids.  For allopolyploids, validate chromosome
assignment and allelic exclusions particularly carefully; separating clear
subgenomes before within-subgenome haplotype phasing may be more appropriate.
