# PHap v2

> [!WARNING]
> This branch is a **testing release**, not a stable production release.
> All 66 automated tests pass, including synthetic triploid, pentaploid, and
> hexaploid dosage paths. On autotetraploid potato, all 48 groups have
> completed HiFi/ONT read assignment, extraction, reassembly, and yak
> evaluation; Hi-C extraction is also complete for all 48 groups. HapHiC
> scaffolding has been validated on two selected groups, but not yet on all 48.
> It is published so that collaborators can evaluate other species and report
> failures, parameter sensitivity, and ploidy-specific behavior. Do not replace
> an established production workflow without validating all QC outputs.

A haplotype-resolved and telomere-to-telomere genome assembly pipeline (PHap) tailored for autopolyploids, relying solely on common sequencing data including long-reads and Hi-C.

PHap v2 adds an auditable allelic-unitig-table workflow. It assigns each unitig
from aggregate chromosome-local PAF evidence, retains rearranged and mixed-
strand segments, uses real reference-coordinate overlap and dosage constraints,
and reports ambiguity explicitly. Parental information is never used for
clustering or read assignment; it may only be used after a run for evaluation.

See [MULTISPECIES_TESTING.md](MULTISPECIES_TESTING.md) before testing another
species. It records the current validation boundary, required metadata, a
small-scale testing strategy, and the files needed for a useful bug report.

See [GENERAL_PLOIDY.md](GENERAL_PLOIDY.md) for dosage labels, depth
classification, non-tetraploid commands, output naming, and required QC.

See [UPGRADE_RECORD.md](UPGRADE_RECORD.md) for the consolidated Chinese record
of code changes, retained and rejected algorithm experiments, current validation
results, and the recommended next run.

See [V2_DEVELOPMENT_VALIDATION.md](V2_DEVELOPMENT_VALIDATION.md) for the
separate development-data, code, parameters, result-directory, and validation
boundary record for the current potato tests.

The complete current-potato validation scripts are
[`examples/run_cluster_full_potato_v2.sh`](examples/run_cluster_full_potato_v2.sh)
and
[`examples/run_phase_reads_full_potato_v2.sh`](examples/run_phase_reads_full_potato_v2.sh).
They use new output directories and stop rather than overwriting an existing
run.

See [ALLELIC_TABLE_V2.md](ALLELIC_TABLE_V2.md) for the algorithm, thresholds,
output interpretation, and validation strategy.

See [CHROMOSOME_EXTRACTION_V2.md](CHROMOSOME_EXTRACTION_V2.md) for the rewritten
`02.chr_seq` assignment rules, consistency checks, and QC outputs.

See [CLUSTER_V2.md](CLUSTER_V2.md) for the rewritten `03.cluster` constraint
solver, dosage handling, Hi-C refinement, and acceptance checks.

See [RECLUSTER_V2.md](RECLUSTER_V2.md) for the rewritten `04.recluster`
confidence rules, iterative propagation, and deferred-sequence handling.

See [RESCUE_V2.md](RESCUE_V2.md) for the rewritten `05.rescue` candidate
scope, chromosome-and-group confidence rules, final partition validation, and
auditable outputs.

See [PHASE_READS_V2.md](PHASE_READS_V2.md) for the staged read-assignment,
collapsed-read balancing, one-pass FASTQ extraction, assembly, and scaffolding
workflow. Variant-based read phasing is a documented later extension.

## PHap v2 quick start

Reuse an existing raw minimap2 PAF:

```shell
python PHap.py allelic_table \
  --p_utg p_utg.fa \
  --mT2T mT2T.fa \
  --contig_type contig_depth.txt \
  --paf p_utg_vs_mT2T.paf \
  --output-dir 02.cluster.v2/01.putg_vs_mT2T
```

Omit `--paf` to run minimap2. The primary result is
`corrected_allelic_table.txt`; review `collinear_chain.qc.tsv`,
`allelic_table.qc.tsv`, `rejected_projections.tsv`, and
`allelic_table.summary.json` before clustering. GFA is optional and is not used
unless `--gfa` is supplied. With `--gfa`, PHap also removes
direct graph-link conflicts and writes `allelic_pairs.gfa_validation.tsv` plus
its JSON summary. The no-sequence hifiasm GFA is sufficient for this check.

Run the complete workflow only through initial clustering with:

```shell
python PHap.py cluster ... --stop_after cluster
```
## Overview
![|600](https://bioin-1320274504.cos.ap-nanjing.myqcloud.com/images/PHAP.overview.v3.png)
1. **Initial assembly.** A primary contig assembly (*p_ctg*) and a phased unitig assembly (*p_utg*) are assembled using hifiasm with the long-read sequencing data including PacBio HiFi and Oxford Nanopore Ultra-Long sequencings. 
2. **mT2T assembly.** A mosaic T2T (mT2T) reference is assembled based on the all-vs-all alignments of the *p_ctg* assembly contigs. 
3. **Allelic unitig table construction.** All unitigs of the *p_utg* assembly are aligned to the mT2T reference to build an allelic table. 
4. **Unitig clustering.** Unitigs in the allelic table are clustered according to the strength of Hi-C signals and guided by alignments against mT2T. 
5. **Unitig re-clustering.** Re-clustering for each of remaining unitigs based on their relative intensity of Hi-C interaction against each group. 
6. **Read phasing and de novo assembly.** Long reads are mapped to the *p_utg* unitigs, assigned into haplotypes, and de novo assembled independently for each haplotype. 
7. **Chromosome-scale assembly.** Contigs of each haplotype assembly are scaffolded with Hi-C data and followed by gap filling and polishing with long reads.

## System Requirements
* All scripts and analyses were developed and tested on Linux operating system environment.
* Several essential tools were required: 
	* [hifiasm](https://github.com/chhylp123/hifiasm)
	* [minimap2](https://github.com/lh3/minimap2)
	* [SAMtools](https://github.com/samtools/samtools)
	* [BWA](https://github.com/lh3/bwa)
	* [SeqKit2](https://github.com/shenwei356/seqkit)
	* [HapHiC](https://github.com/zengxiaofei/HapHiC)
	* [PanDepth](https://github.com/HuiyangYu/PanDepth)
	* [Mash](https://github.com/marbl/Mash?tab=readme-ov-file)
	* [TGS-GapCloser](https://github.com/BGI-Qingdao/TGS-GapCloser)
	* [Winnowmap2](https://github.com/marbl/Winnowmap)
	* [T2T-polish](https://github.com/arangrhie/T2T-Polish "")
	* [Python 3.9.7](https://www.python.org/downloads/)

## Installation
* Download the current multispecies testing branch
```shell
$ git clone --branch testing/v2-multispecies-20260823 \
    https://github.com/pxxiao-hz/PHap.v2.git
$ cd /path/to/PHap.v2/
$ python -m pip install -r requirements.txt
$ chmod +x PHap.py
```

## Usage
* View the help document of PHap
```shell
$ /path/to/PHap/PHap.py -h

  Usage: phap [command] <parameters>

  Command       Description
  --------      ------------------------------------------------------------------------
  mt2t          Generate a mosaic telomere-to-telomere genome using the p_ctg genome as
                the reference. This is used for creating an allelic contig table to solve
                the allelic conflict problem.

  cluster       Cluster contigs. First, extract contigs corresponding to each chromosome
                based on mT2T alignment. Then, cluster the contigs for each haplotype of
                each chromosome using Hi-C signals. Third, cluster the contigs not on mT2T
                chromosomes to clustered group.

  phase_reads   Based on the clustering results, the corresponding reads of each haplotype
                are extracted, then assembled and anchored separately.

  Use phap [command] --help/-h to see detailed help for each individual command.
```
* View the help document of subcommands of PHap
```
$ /path/to/PHap/PHap.py mt2t -h 
usage: Get mosaic T2T (mT2T) reference from primary contig assembly (p_ctg).

$ /path/to/PHap/PHap.py cluster -h
usage: Haplotype clustering of autopolyploid genome.

$ /path/to/PHap/PHap.py phase_reads -h
usage: Haplotype assembly and scaffolding of autopolyploid genome.
```
* For convenience, users can add `/path/to/PHap/PHap.py` to environment variables, such as: `ln -s /path/to/PHap/PHap.py ~/.local/bin/phap`, and then the `phap`command can be executed directly from any location in the terminal.

## The pipeline for assembling a tetraploid potato genome
Please check the [Pipeline](Pipeline.md).

The upgraded read-assignment, reassembly, and scaffolding design is documented
in [PHASE_READS_V2.md](PHASE_READS_V2.md).

If HiFi/ONT assembly and Hi-C extraction were completed in separate runs, use
`phap phase_reads --scaffold-only --assembly-dir DIR --hic-reads-dir DIR` to
run HapHiC directly from those existing results. This mode supports
`--groups`/`--chromosomes`, `--temp-dir`, and per-group `--resume` checkpoints;
it does not require BAM inputs or `--contig-type`.

## Note

PHap now implements integer dosage and dynamic group counts for ploidies of at
least two. Synthetic tests cover non-tetraploid paths through the v2 workflow.
Real-data development and end-to-end calibration remain focused on
autotetraploid potato (_Solanum tuberosum_), so another species or ploidy must
still be treated as experimental and validated as described in
[GENERAL_PLOIDY.md](GENERAL_PLOIDY.md) and
[MULTISPECIES_TESTING.md](MULTISPECIES_TESTING.md).
