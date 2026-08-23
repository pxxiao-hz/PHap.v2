# Validation on the Current Workspace Data

Run date: 2026-08-20

Inputs:

- `hs.100k.asm.bp.p_utg.gfa.remove.plastid.contamination.fa`
- `hs.mT2T.v2.fa`
- `contig_depth.txt`
- `02.cluster/01.putg_vs_mT2T/p_utg_vs_mT2T.paf`
- `hs.100k.asm.bp.p_utg.noseq.gfa`

Command:

```shell
python PHap.v2/PHap.py allelic_table \
  --p_utg hs.100k.asm.bp.p_utg.gfa.remove.plastid.contamination.fa \
  --mT2T hs.mT2T.v2.fa \
  --contig_type contig_depth.txt \
  --paf 02.cluster/01.putg_vs_mT2T/p_utg_vs_mT2T.paf \
  --gfa hs.100k.asm.bp.p_utg.noseq.gfa \
  --output-dir 02.cluster.v2/01.putg_vs_mT2T
```

Results:

| Metric | Value |
| --- | ---: |
| Raw query unitigs in PAF | 4,416 |
| Accepted collinear chains | 3,094 |
| Dense accepted chains | 3,067 |
| Sparse-syntenic accepted chains | 27 |
| Final localized projections | 2,664 |
| Final table rows | 1,754 |
| Final table unitigs | 1,548 |
| Allelic pairs | 1,196 |
| Pairs with non-positive coordinate overlap | 0 |
| Rows with dosage greater than four | 0 |
| GFA direct-link pairs remaining | 0 |
| GFA bubble-like supported pairs | 202 |
| GFA shared-read review pairs | 149 |
| GFA-neutral pairs | 845 |
| GFA-conflict atomic intervals omitted | 1,344 (25,182,403 bp) |
| Ambiguous over-capacity intervals omitted | 30 |

`utg000006l` and `utg000018l` demonstrate why aligned-base coverage and anchor
span must be separated. Their aligned-base coverages are only 24.28% and
23.90%, but their high-MAPQ collinear anchors span 98.03% and nearly 100% of
their query lengths. Both have about 95.1% chain identity and unambiguous chr11
placement. They are retained as `sparse_syntenic`; the final table contains
47.84 Mb and 30.02 Mb for them, with 30.02 Mb shared.

The GFA is neutral for `utg000006l`/`utg000018l`: their graph degrees are zero
and two, respectively, but they have no direct edge, common graph neighbor, or
shared `A`-record read. Thus the GFA does not reject this pair, while its
allelic assignment remains supported mainly by the long, unambiguous chr11
collinear projections rather than by graph topology.

Before graph-aware filtering, 1,234 of 2,728 proposed allelic pairs were direct
GFA neighbors. Every one of those pairs overlapped by less than 50 kb in the
table, and the table overlap was bounded by the GFA link overlap (allowing a
20 kb projection tolerance). This is strong evidence that adjacent assembly
segments had been misinterpreted as alleles. The graph-aware rerun removes all
direct-link pairs from the final table. It does not delete any unitig sequence
from the assembly FASTA.

In the original corrected 100 kb table, `utg000006l` occurred in 476 bins and
`utg000018l` in 301 bins. All 301 `utg000018l` bins also contained
`utg000006l`, and every shared bin contained four unitigs. Across that 30.1 Mb
shared interval, the five other unitigs were `utg000086l`, `utg000196l`,
`utg000281l`, `utg000457l`, and `utg001302l`.

These are internal consistency checks, not a biological precision estimate.
The next validation layer should use parental markers, known phased blocks, or
simulated truth to estimate precision and recall.

## Chromosome Sequence Extraction

The rewritten extractor was run against the accepted collinear PAF and its QC:

```shell
python PHap.v2/utils/extract_chr_from_putg.py \
  --p_utg hs.100k.asm.bp.p_utg.gfa.remove.plastid.contamination.fa \
  --mT2T hs.mT2T.v2.fa \
  --paf 02.cluster.v2/01.putg_vs_mT2T/putg_vs_mT2T.collinear.paf \
  --chain-qc 02.cluster.v2/01.putg_vs_mT2T/collinear_chain.qc.tsv \
  --wd 02.cluster.v2/02.chr_seq --chr_num 12 \
  --min-chromosome-margin 0.05
```

| Metric | Value |
| --- | ---: |
| Input unitigs | 4,564 |
| Input sequence | 2,580,736,761 bp |
| Strict-chain assigned unitigs | 3,094 |
| Best-chromosome rescued unitigs | 1,139 |
| Total assigned unitigs | 4,233 (2,553,586,569 bp) |
| Unassigned unitigs | 331 (27,150,192 bp) |
| Assigned fraction of p_utg sequence | 98.95% |
| Output unique unitigs | 4,564 |
| Output sequence | 2,580,736,761 bp |
| Duplicate output IDs | 0 |

The input and output record and base counts are identical. `utg000006l` and
`utg000018l` are both assigned to chr11 as `sparse_syntenic`, with query anchor
spans of 98.03% and 100%, respectively.

Low allelic-projection coverage no longer removes sequence from a chromosome.
For example, `utg000010l`, `utg000014l`, and `utg000017l` are rescued to chr09,
chr12, and chr07 because their best-target margins are 89.26%, 74.03%, and
91.90%, respectively.

The 331 unassigned unitigs comprise 148 with no alignment record, 137 with no
primary alignment after the basic alignment filters, and 46 with no reliable
chromosome target or a best-target margin below 5%. Only one unassigned unitig
is longer than 1 Mb: `utg001680l` (1,587,195 bp), which has no accepted
collinear-chain record.

## Constraint-Aware Clustering

The rewritten clusterer was run independently for chr01 through chr12 with
ploidy four, `balance_weight=1.0`, ten refinement rounds, and default `weakest`
constraint relaxation.

| Metric | Value |
| --- | ---: |
| Allelic-table unitigs assigned | 1,548 / 1,548 |
| Dosage group memberships | 2,381 |
| Dosage/group-count mismatches | 0 |
| Enforced allelic conflicts | 0 |
| Unassigned table unitigs | 0 |
| Chromosomes with a failed constraint search | 0 |
| Relaxed allelic constraints | 1 |

The sole relaxation is `utg000079l`/`utg000558l` on chr01. These unitigs occur
in a five-haplotig clique whose summed dosage is five at ploidy four, making the
unmodified conflict graph mathematically infeasible. This edge has only 12,038
bp of table overlap and a minimum directional overlap ratio of 2.29%, so it is
the weakest edge in that clique. The decision and final shared group are stored
in `chr01/cluster_relaxed_constraints.tsv`.

The four chr11 seed groups contain 47.38, 52.60, 44.22, and 14.01 Mb. This
imbalance should not be corrected by deleting `utg000006l` or `utg000018l`:
four mutually exclusive long unitigs of 47.34, 34.12, 29.35, and 6.24 Mb anchor
the four groups. It reflects the conservative seed content and haplotype length
difference; non-table chromosome unitigs are handled in `04.recluster`.

## Chromosome Unitig Reassignment

The rewritten `04.recluster` was run on all twelve chromosome FASTA files with
dosage-adjusted Hi-C, a minimum adjusted support of 5, a density margin of
0.10, and low-confidence policy `defer`.

| Metric | Value |
| --- | ---: |
| Chromosome input unitigs | 4,233 (2,553,586,569 bp) |
| Fixed allelic-table seeds | 1,548 |
| High-confidence Hi-C additions | 889 |
| Dosage-four additions | 29 |
| Total assigned | 2,466 (2,469,506,203 bp; 96.71%) |
| Explicitly deferred | 1,767 (84,080,366 bp; 3.29%) |
| No assigned-neighbor Hi-C | 1,361 |
| Insufficient adjusted links | 340 |
| Group-margin below 0.10 | 66 |
| Assigned dosage errors | 0 |
| Changed seed assignments | 0 |
| Partition errors | 0 |

The dosage table lacks records for 885 chromosome unitigs totaling 36,023,057
bp. They are marked `default_haplotig`; 780 of them (29,341,834 bp) remain
deferred. A strict run can use `--recluster_unknown_dosage_policy error`.

For chr11, the final groups contain 47.81, 52.64, 50.69, and 52.15 Mb.
`utg000281l`, `utg000196l`, and `utg000005l` join g4, while `utg000075l`
joins g3. All four decisions have strong normalized Hi-C margins.

The 6.5 GB CLM contained 2,416,352 records. The new splitter scanned it once in
11.6 seconds with about 40 MB peak memory and generated all 48 nonempty group
CLM files. The old design would have scanned the complete CLM independently
for every chromosome.

## Chromosome-Unassigned Unitig Rescue

The rewritten `05.rescue` was run with minimum adjusted support 5, chromosome
and dosage-boundary group margins 0.10, unsupported dosage policy `defer`, and
low-confidence policy `defer`.

| Metric | Value |
| --- | ---: |
| Original p_utg input | 4,564 unitigs; 2,580,736,761 bp |
| Immutable step 04 assignments | 2,466 unitigs; 2,469,506,203 bp |
| Immutable step 04 deferrals | 1,767 unitigs; 84,080,366 bp |
| Actual rescue candidates (`un_chr.fa`) | 331 unitigs; 27,150,192 bp |
| High-confidence rescued | 41 unitigs; 6,140,819 bp |
| Rescue candidates still unassigned | 290 unitigs; 21,009,373 bp |
| Final assigned | 2,507 unitigs; 2,475,647,022 bp (95.93%) |
| Final explicitly unassigned | 2,057 unitigs; 105,089,739 bp (4.07%) |
| Final dosage memberships | 3,622 |
| Dosage, changed-step-04, partition, or length errors | 0 |

The 41 rescued unitigs comprise 35 haplotigs, five diplotigs, and one
triplotig. All meet the configured thresholds; the observed minima are 5.0
adjusted links, 0.202 chromosome margin, and 0.163 dosage-boundary group
margin.

The 290 deferred candidates comprise 56 with no assigned-neighbor Hi-C, 34
with insufficient adjusted links, nine with group margin below 0.10, 166
`other` unitigs, and 25 `replotig` unitigs. Unsupported types are not guessed
as haplotigs by default.

The final FASTA audit found 2,507 unique assigned IDs with exactly 3,622
dosage memberships and 2,057 unique unassigned IDs. The two sets are disjoint
and their union is all 4,564 p_utg unitigs; every emitted sequence length
matches the original assembly.

The rescue run completed in 37.7 seconds with about 317 MB peak memory. The
6.5 GB CLM was then scanned once in 13.6 seconds with about 44 MB peak memory,
writing 235,948 final group-record memberships across 48 CLM files.
