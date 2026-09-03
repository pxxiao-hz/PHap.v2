# PHap v2 Read Phasing, Assembly, and Scaffolding Upgrade

Updated: 2026-08-31

## Scope

This upgrade covers read assignment, FASTQ extraction, per-group hifiasm
assembly, and HapHiC scaffolding after `05.rescue`. It must not use parental
information. Parental k-mers may be used only after a run to evaluate the
method.

Variant-based read phasing is deliberately deferred. The first implementation
uses alignment evidence where it is informative and an auditable balanced
fallback where reads from a collapsed unitig are biologically
indistinguishable.

Implementation status: stages 1–6 are implemented in
`utils/phase_reads_assignment.py` and
`utils/phase_reads_assemble_anchor.py`. Synthetic validation covers dosage
checking, deterministic collapsed-read balancing, pair-level Hi-C decisions,
single-pass FASTQ dispatch, BAM-offset and data-type checkpoints, per-group
assembly/scaffolding checkpoints, bounded parallel assembly, output validation,
and preservation of an earlier assembly when hifiasm fails. Coordinate-sorted
BAM evidence for HiFi/ONT is aggregated through SQLite instead of retaining
all long reads in memory. Hi-C now defaults to a corrected v1 streaming path:
primary read1 plus RNEXT supplies the two mate references, collapsed pairs are
assigned mutually exclusively, and no Hi-C evidence SQLite is created. The
older detailed Hi-C pair-evidence implementation remains available through
`--hic-assignment-backend sqlite`.

All 66 automated tests pass, including six-group dosage validation. Real
autotetraploid potato data have completed
HiFi/ONT assignment, extraction, and hifiasm reassembly for all 48 groups;
all 48 assemblies have also been evaluated with yak. Corrected fast-v1 Hi-C
assignment and paired FASTQ extraction completed for all 48 groups. The
scaffold-only interface has completed on two selected groups
(`chr03_group1` and `chr09_group2`), including BWA mapping, BAM filtering,
HapHiC clustering/reassignment/sorting/building, and per-group checkpoints.
Full 48-group scaffolding and Juicebox output have deliberately not been run,
so this remains a testing release rather than a complete production
validation.

Use `--temp-dir DIR` to place PHap working databases/directories, Python
temporary files, SQLite temporary files, and temporary files created by child
tools under a chosen location. Final outputs still go to `--output-dir`.
Completed temporary evidence databases and staging directories are removed;
interrupted evidence databases are retained there for checkpoint resume.

## Staged implementation

### Stage 1: validated group and dosage model

- Parse the final `group.reassignment.cluster.txt` in both directions:
  `group -> unitigs` and `unitig -> groups`.
- Treat final group membership as the effective dosage of an assigned unitig.
- Cross-check supported contig types without silently converting unsupported
  values to haplotigs.
- Require a complete, contiguous set of haplotype groups for every chromosome.

### Stage 2: read-centric assignment

- Filter unmapped, secondary, duplicate, QC-failed, short, and low-identity
  alignments. Retain primary alignments at every MAPQ because highly similar
  polyploid haplotypes can produce accurate but non-unique placements.
- Record MAPQ for auditing, but do not use it to filter or weight phase-read
  evidence. Alignment scores use query coverage and identity.
- Retain useful supplementary alignments and union overlapping query evidence.
- Aggregate alignment evidence by chromosome-haplotype group.
- Assign a read to one group only when the best-vs-second evidence margin
  passes the configured threshold.
- Defer contradictory or weak reads explicitly.
- Stream BAM records into a resumable SQLite evidence store, aggregate one read
  at a time from disk, write the final decision index, and release each
  technology before starting the next one.

### Stage 3: balanced collapsed-read fallback

- Apply the fallback only when a read has no group-unique alignment evidence
  and its strongest groups belong to one chromosome.
- Partition reads exclusively across their candidate groups: never give the
  same read to every group and never silently discard the unused fraction.
- Sort reads deterministically and greedily assign them to the candidate group
  with the lowest estimated read depth.
- Estimate depth from assigned read bases divided by the total unitig bp in
  each group, so the fallback supports comparable assembly coverage rather
  than merely equal read counts.
- Record every fallback assignment and the final per-group depth estimate.

### Stage 4: Hi-C pair assignment

- The default `--hic-assignment-backend fast-v1` scans only primary read1
  records and uses BAM RNEXT for the mate reference. It filters no record by
  MAPQ and creates no evidence or assignment SQLite.
- Intersect the existing group memberships of the two complete unitigs. Accept
  a unique intersection and defer incompatible mates.
- Balance collapsed-only candidates mutually exclusively among only the groups
  already assigned to their complete unitigs; never create local regions or
  blocks.
- Write one read-name file per group for seqkit extraction. The detailed legacy
  v2 scorer is retained as `--hic-assignment-backend sqlite` for controlled
  comparison, not as the large-data default.

### Stage 5: one-pass FASTQ dispatch

- Scan HiFi FASTQ once and ONT FASTQ once.
- In the detailed SQLite mode, scan Hi-C R1/R2 together once and validate mate
  names. In fast-v1 mode, use each group's read-name list with seqkit; this
  reproduces the v1 extraction architecture and validates R1/R2 output counts.
- Dispatch each accepted record directly to its group output.
- Report requested, found, missing, duplicated, and unassigned records.
- For HiFi and ONT, `--extract-backend auto` selects the high-throughput
  `seqkit fx2tab -> gawk assignment demultiplexing -> pigz` backend when all
  three tools are available. The assignment map is exported sequentially from
  SQLite, the raw FASTQ is scanned once, and named pipes feed one single-thread
  pigz stream per output group without writing uncompressed intermediate
  FASTQ files. fast-v1 Hi-C extraction uses `seqkit grep -> gawk count -> pigz`
  with bounded group concurrency and writes no uncompressed FASTQ intermediate.
- The fast backend preserves the full read header, sequence, and quality
  string. FASTQ plus lines are normalized to `+`. When an ONT mean-quality
  threshold is requested, gawk computes the same arithmetic Phred mean as the
  Python implementation; the calculation is skipped entirely at the default
  threshold of zero.

### Stage 6: reliable assembly and scaffolding

- Separate `assign`, `extract`, `assemble`, and `scaffold` stages.
- Use bounded `jobs x threads_per_job` resource scheduling.
- Propagate every external command failure and validate expected outputs before
  starting the next stage.
- Resolve tools from `PATH` or explicit CLI paths; remove user-specific paths.
- Write an input/parameter manifest and stage summaries for safe resumption.
- Checkpoint HiFi, ONT, and Hi-C assignment/extraction independently and
  checkpoint every assembly/scaffolding group independently.
- Log stage progress, record counts/rates, and group completion percentages to
  stderr and `phase_reads.log`.
- `--resume --rerun-from STAGE` keeps checkpoints before `STAGE` and actively
  reruns that stage plus downstream stages. Changing jobs/threads does not
  invalidate biological-result checkpoints.

### Small-scope validation

- `--data-types hifi` (or `ont`/`hic`) permits assignment testing as individual
  BAM files become available. Inputs are required only for selected data types;
  full assembly still requires HiFi and scaffolding requires Hi-C.
- `--chromosomes chr01 chr02` runs all haplotype groups for selected
  chromosomes.
- `--groups chr01_group1 chr02_group3` writes only those groups. Assignment
  still includes every sibling haplotype group on the corresponding
  chromosomes so collapsed-read balancing is not biased toward the requested
  outputs.
- Comma-separated values are also accepted. `--groups` and `--chromosomes` are
  mutually exclusive.

### Scaffold-only mode

When assembly and Hi-C extraction were completed in separate runs, use
`--scaffold-only` to connect the two immutable result directories without
rerunning assignment, extraction, or assembly:

```shell
phap phase_reads \
  --scaffold-only \
  --group 02.cluster.v2/05.rescue/group.reassignment.cluster.txt \
  --assembly-dir 03.phase_reads.hifi_ont/03.assembly \
  --hic-reads-dir 03.phase_reads.hic/02.reads \
  --output-dir 03.phase_reads.scaffold \
  --temp-dir 03.phase_reads.scaffold.tmp \
  --groups chr03_group1 chr09_group2 \
  --jobs 2 \
  --threads-per-job 10 \
  --haphic-processes 1
```

The assembly directory must contain
`GROUP/GROUP.asm.bp.p_ctg.gfa`; the Hi-C directory must contain
`GROUP.Hi-C.1.fq.gz` and `GROUP.Hi-C.2.fq.gz`. Final results are written under
`OUTPUT/04.scaffold`. The external input files are opened read-only and are
recorded by path, size, and modification time in `run_manifest.json`.

Scaffold-only mode does not require BAM files, raw whole-dataset FASTQs, or
`--contig-type`. It supports `--groups`, `--chromosomes`, `--temp-dir`, and
per-group scaffold checkpoints. Repeat the same command with `--resume` to
skip completed groups. If a selected GFA or Hi-C FASTQ changes, resume is
rejected instead of silently reusing a stale scaffold. Use
`--resume --rerun-from scaffold` to deliberately rebuild the selected groups.

### Deferred stage: variant phasing

Later, heterozygous variants and read phase blocks can replace balanced
fallback assignments when phase blocks can be connected to group-specific
unitigs. This is not part of the current upgrade.

## Acceptance criteria

1. A HiFi or ONT read occurs in at most one group unless a future explicit
   sharing policy is selected.
2. Collapsed-only reads assigned by the fallback are conserved, mutually
   exclusive, deterministic, and balanced by estimated depth.
3. Hi-C mates are always retained or deferred together.
4. FASTQ extraction counts agree with the assignment tables.
5. Changing an input or assignment parameter cannot silently reuse stale
   cached results.
6. A failed extraction, assembly, mapping, or scaffolding command fails the
   corresponding stage and prevents downstream execution.
7. HiFi/ONT and detailed-mode Hi-C FASTQ are scanned once per extraction run;
   fast-v1 Hi-C intentionally uses the v1 group-pattern extraction strategy.
8. All decisions are reproducible without parental information.
9. An interrupted BAM scan resumes from its committed virtual offset; an
   interrupted downstream run restarts only the unfinished technology or
   group.
10. The fast extraction backend must produce the same group, read header,
    sequence, and quality assignments as the Python backend and must fail on
    missing, duplicated, malformed, or compressor-error conditions.

### Extraction backend controls

- `--hic-assignment-backend fast-v1|sqlite` defaults to `fast-v1`. `fast-v1`
  requires `seqkit`, `gawk`, and `pigz` when extraction is requested.

- `--extract-backend auto|seqkit|python` defaults to `auto`. `auto` uses the
  fast backend when `seqkit`, `gawk`, and `pigz` are on `PATH`, and otherwise
  logs a warning and uses Python.
- `--extract-threads 8` controls parallel gzip input decoding in seqkit. Each
  group output is compressed by a separate one-thread pigz process.
- `--seqkit`, `--gawk`, and `--pigz` accept explicit executable paths.
- Backend and thread choices are recorded as runtime manifest information and
  do not invalidate an already completed biological extraction checkpoint.

## Remaining real-data validation

For additional species, start with clear chromosomes and difficult/repetitive
chromosomes. Compare read conservation, fallback fraction, group depth,
assembly size/N50, read mapping, Hi-C valid-pair rate, scaffold contact maps,
runtime, memory, and I/O. For the current potato data, the remaining optional
validation is full 48-group scaffolding and Juicebox/contact-map review.
Parental yak evaluation is performed only after all method decisions have been
fixed.

The current alignment-length, identity, group-margin, and balancing-tolerance
defaults are safe initial values rather than final biological thresholds. They
must be calibrated from assignment summaries and mappings on the supplied real
data without using parental labels for parameter selection. MAPQ is retained
only as an audit field and is not a phase-read filter or score weight.
