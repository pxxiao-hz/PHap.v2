# PHap v2 multispecies testing guide

## Release status

Branch: `testing/v2-multispecies-20260823`

This is a validation branch. All 58 automated tests pass. On real
autotetraploid potato data, all 48 groups have completed HiFi/ONT assignment,
extraction, reassembly, and yak evaluation, and all 48 paired Hi-C group FASTQs
have been extracted. HapHiC scaffolding has completed on two selected groups;
full 48-group scaffolding and Juicebox review were not run.

Parental information must not be supplied to clustering or read assignment.
Parental k-mers or markers may be used only after the workflow is fixed, to
evaluate accuracy.

## Important current limitations

- Development and current real-data calibration focus on autotetraploid
  potato.
- `--top_n` specifies expected haplotype count and `--chr_num` specifies the
  chromosome count, but general triploid, hexaploid, and higher-ploidy support
  is not yet certified.
- Thresholds were not calibrated on diverse genome sizes, repeat content,
  heterozygosity levels, or sequencing depths.
- `phase_reads` has extensive potato validation, but only two groups have been
  scaffolded and no other species has completed end-to-end validation. Treat
  results on a new species as experimental.
- A chromosome- or group-limited `phase_reads` run still scans the complete raw
  FASTQ because FASTQ is not chromosome-indexed.

## Recommended test sequence

1. Record species, expected ploidy, chromosome count, estimated genome size,
   sequencing technologies, and approximate depths.
2. Run the bundled tests before using biological data:

   ```shell
   python -m unittest discover -s tests -v
   ```

3. Start with allelic-table generation and inspect every QC TSV/JSON before
   clustering.
4. Run one chromosome through `cluster` using the correct `--top_n` and
   `--chr_num`. Use a new output directory; do not overwrite a previous run.
5. Review `03.cluster`, `04.recluster`, and `05.rescue` validation summaries,
   deferred unitigs, relaxed constraints, group lengths, and Hi-C support.
6. Test `phase_reads` on one chromosome with `--chromosomes`, one sequencing
   type with `--data-types`, and a dedicated `--temp-dir`.
7. Expand to all chromosomes only after small-scale QC is acceptable.

Example phase_reads engineering test:

```shell
python PHap.py phase_reads \
  --bam-hifi /path/to/hifi_to_putg.sort.bam \
  --data-types hifi \
  --contig-type /path/to/contig_depth.txt \
  --group /path/to/group.reassignment.cluster.txt \
  --chromosomes chr01 \
  --output-dir phase_reads.chr01 \
  --temp-dir phase_reads.tmp \
  --stop-after assign
```

After inspecting assignment QC, add the raw HiFi FASTQ and resume extraction:

```shell
python PHap.py phase_reads \
  --bam-hifi /path/to/hifi_to_putg.sort.bam \
  --hifi /path/to/hifi.fastq.gz \
  --data-types hifi \
  --contig-type /path/to/contig_depth.txt \
  --group /path/to/group.reassignment.cluster.txt \
  --chromosomes chr01 \
  --output-dir phase_reads.chr01 \
  --temp-dir phase_reads.tmp \
  --resume --stop-after extract
```

Input BAM, FASTQ, group, and contig-type files are opened read-only. Temporary
work is written under `--temp-dir`; final results are written under
`--output-dir`.

## What to report

Please include:

- the exact Git commit and complete command line;
- species, expected ploidy/chromosome count, genome size, and data depths;
- Python and external-tool versions;
- `phase_reads.log` or the relevant stage log;
- stage summary JSON files and small QC TSV files;
- peak memory, elapsed time, and the stage where a failure occurred;
- whether the failure is reproducible from a clean output directory.

Do not upload raw reads, BAM files, assemblies, credentials, private paths, or
other sensitive data to a public issue unless sharing them is explicitly
authorized.
