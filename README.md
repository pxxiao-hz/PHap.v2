# PHap
A haplotype-resolved and telomere-to-telomere genome assembly pipeline (PHap) tailored for autopolyploids, relying solely on common sequencing data including long-reads and Hi-C.
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
	* [Python 3.9+](https://www.python.org/downloads/)

External executables must be available on `PATH`. PHap checks required tools
before starting each expensive workflow stage. HapHiC's `haphic` and
`filter_bam` commands must also be exposed on `PATH`.

## Installation

```shell
git clone https://github.com/pxxiao-hz/PHap.v2.git
cd PHap.v2
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
phap --help
```

For development, install `python -m pip install -e ".[dev]"`.

## Usage

```shell
phap --help
phap --version
```

View subcommand-specific options:

```
phap mt2t --help
usage: Get mosaic T2T (mT2T) reference from primary contig assembly (p_ctg).

phap cluster --help
usage: Haplotype clustering of autopolyploid genome.

phap phase_reads --help
usage: Haplotype assembly and scaffolding of autopolyploid genome.
```

`python PHap.py ...` remains available as a source-checkout compatibility
entry point.

## The pipeline for assembling a tetraploid potato genome
Please check the [Pipeline](Pipeline.md).

## Note
**PHap** is currently developed to do _haplotype-resolved telomere-to-telomere (T2T) genome assembly_ in **autotetraploid genomes**, with a primary focus on crops such as potato (_Solanum tuberosum_). 
Theoretically, **PHap can be easily extended to support other polyploid genome types**, including **triploid**, **hexaploid**, and more complex genomes.
