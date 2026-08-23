# PHap v2 Read Phasing, Assembly, and Scaffolding Upgrade

Updated: 2026-08-22

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
BAM evidence is aggregated through SQLite instead of retaining all reads and
all three sequencing technologies in memory. A chr01 HiFi assignment run has
completed and validated BAM-offset resume and depth balancing. HiFi FASTQ
extraction is currently being tested; ONT/Hi-C extraction, reassembly, and
scaffolding have not yet completed real-data validation.

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

- Filter unmapped, secondary, duplicate, QC-failed, low-MAPQ, short, and
  low-identity alignments.
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

- Score R1 and R2 separately, then make one pair-level decision.
- Accept pairs whose ends agree, or whose unique end is compatible with the
  ambiguous end.
- Defer cross-chromosome or cross-haplotype conflicts.
- Balance collapsed-only ambiguous pairs without separating their mates.

### Stage 5: one-pass FASTQ dispatch

- Scan HiFi FASTQ once and ONT FASTQ once.
- Scan Hi-C R1/R2 together once and validate mate names.
- Dispatch each accepted record directly to its group output.
- Report requested, found, missing, duplicated, and unassigned records.

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
7. Each large FASTQ input is scanned once per extraction run.
8. All decisions are reproducible without parental information.
9. An interrupted BAM scan resumes from its committed virtual offset; an
   interrupted downstream run restarts only the unfinished technology or
   group.

## Real-data validation after data are supplied

Start with clear chromosomes such as chr01/chr11 and difficult chromosomes
chr02/chr05/chr08. Compare read conservation, fallback fraction, group depth,
assembly size/N50, read mapping, Hi-C valid-pair rate, scaffold contact maps,
runtime, memory, and I/O. Parental yak evaluation is performed only after all
method decisions have been fixed.

The current MAPQ, alignment-length, identity, group-margin, and balancing
tolerance defaults are safe initial values rather than final biological
thresholds. They must be calibrated from assignment summaries and mappings on
the supplied real data without using parental labels for parameter selection.
