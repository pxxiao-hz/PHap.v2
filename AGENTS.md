# PHap Development Constraints

This file applies to the entire repository. Its purpose is to keep upgrades
reproducible, scientifically auditable, and safe for large genome-assembly
workloads.

## 1. Upgrade priorities

Work in this order unless a task explicitly requires otherwise:

1. Make the command-line interface installable and runnable on a clean system.
2. Add deterministic tests around the current scientific behavior.
3. Fix correctness and data-loss defects.
4. Remove duplication and split orchestration from algorithms.
5. Optimize runtime and memory only after correctness is measured.
6. Change scientific heuristics only with before/after assembly evidence.

Do not combine a behavior-preserving refactor with a scientific algorithm
change in the same commit.

## 2. Scientific invariants

- Never silently drop, duplicate, rename, reverse-complement, or merge a
  contig, unitig, or read. Every such operation must be represented in an
  audit table with the source ID, destination ID/group, reason, and evidence.
- Preserve FASTA/FASTQ identifiers exactly unless an output format explicitly
  requires a new ID. When IDs are changed, emit a two-column mapping file.
- Treat PAF coordinates as 0-based half-open intervals. Keep BAM coordinates
  in their native convention internally and document any conversion at the
  boundary. Do not mix coordinate systems in one data structure.
- Parse PAF optional tags by tag name (for example, `tp:A:P`), never by a
  fixed optional-field position.
- Coverage and match ratios must use interval-union coverage when alignments
  can overlap. Summing overlapping alignment lengths is not valid coverage.
- Count each Hi-C pair once. State and test the filtering policy for unmapped,
  secondary, supplementary, duplicate, low-MAPQ, and non-proper alignments.
- The ploidy/haplotype count is one explicit parameter used throughout the
  pipeline. Do not hard-code four groups, tetraploid dosage limits, 12
  chromosomes, or chromosome names in reusable logic.
- Dosage-aware assignment must conserve copy count. A unitig with dosage `d`
  may occupy exactly `d` haplotype groups unless it is explicitly reported as
  ambiguous or unassigned.
- Collapsed-unitig read partitioning must be deterministic and disjoint across
  destination haplotypes. A fixed seed must reproduce the same complete
  partition; it must not select the same subset independently for every group.
- Ambiguous or unsupported records go to explicit `ambiguous`/`unassigned`
  outputs. Lack of Hi-C links is not evidence that a sequence is false.
- The same inputs, parameters, tool versions, and seed must produce
  byte-stable tabular outputs. Sort sets and dictionary-derived output with
  documented tie-breakers.

## 3. Architecture and code boundaries

- Keep `PHap.py` as a thin CLI dispatcher. New scientific logic belongs in
  importable modules, not command wrappers.
- Separate pure parsing/scoring/assignment functions from filesystem changes
  and external process execution.
- Shared FASTA, PAF, BAM, grouping, logging, and subprocess helpers must have
  one canonical implementation. Do not copy helpers between scripts.
- Do not add another monolithic script. Prefer small modules with typed inputs
  and explicit return values.
- Avoid wildcard imports, mutable process-wide state, and `os.chdir()`.
  Resolve paths with `pathlib.Path` and pass working directories explicitly.
- Library functions must not print per-record debug output. Use structured
  logging; default CLI output should remain concise.
- Delete obsolete implementations and commented-out code after replacement is
  covered by tests. Git history is the archive.

## 4. CLI and workflow contract

- Every documented subcommand and `--help` must work after installation
  without importing optional heavy dependencies unnecessarily.
- A user-provided option must either affect execution or be rejected. Never
  expose an option and then replace it with a hard-coded value.
- Validate all input paths, numeric ranges, FASTA IDs, required columns, and
  external executables before starting expensive work.
- Support an explicit output directory. Do not write fixed filenames into the
  caller's current directory without documenting them.
- Return non-zero on any failed stage, including failures in worker processes.
  Do not log an error and continue with partial output.
- Do not use `nohup` inside the application. Job detachment belongs to the
  workflow scheduler or the user.
- Resuming based only on file existence is forbidden. A resumable stage needs
  an atomic completion marker/manifest containing input fingerprints,
  parameters, PHap version, and external tool versions.
- Write outputs to temporary files and atomically rename them after successful
  completion. Never treat a partial file as a valid cache.

## 5. External tools and resources

- Do not add personal absolute paths, home-directory tool paths, or conda
  activation commands. Discover tools from `PATH` or accept explicit CLI/config
  paths.
- Prefer `subprocess.run([...], check=True)` with an argument list. Open
  redirection files in Python. Use a shell only for a documented,
  unavoidable pipeline and quote every user-controlled value.
- Record the exact command and tool version for each stage without recording
  secrets or full read data.
- Bound concurrency with one resource model. `processes * threads_per_process`
  must not exceed the configured CPU budget unless the user opts in.
- Close pools, BAM files, and output handles with context managers. Worker
  exceptions must be collected and re-raised.

## 6. Assembly and clustering changes

- mT2T joins must be based on a conflict-aware oriented overlap graph, not
  independent pairwise concatenation. A contig may be consumed by at most one
  accepted path; cycles, branches, and conflicting orientation must be
  reported rather than guessed.
- Sequence joins must derive trim points from validated alignment chains and
  preserve unconsumed sequence. Test forward and reverse joins, contained
  contigs, repeats, single-alignment cases, and competing overlaps.
- Clustering scores must state how Hi-C counts are normalized for restriction
  sites, contig length, dosage, and group size. Returned normalized values must
  actually be used.
- Assignment thresholds need CLI/config exposure, units, defaults, and a
  sensitivity test. Dataset-specific contig IDs and debug thresholds are not
  allowed in production logic.
- Read phasing must resolve multi-mapping explicitly and report per-group read
  counts, cross-group overlaps, unassigned reads, and estimated contamination.
- Changes intended to improve assembly quality must compare at least
  completeness, continuity, switch/misjoin evidence, haplotype duplication,
  and read/Hi-C support. Do not optimize N50 alone.

## 7. Tests and verification

- Add `pyproject.toml` with pinned/compatible Python dependencies and
  development tools before the first functional refactor.
- Put tests under `tests/`; use tiny synthetic FASTA, PAF, grouping, and SAM/BAM
  fixtures that can run without production datasets.
- Every bug fix needs a regression test that fails on the old behavior.
- At minimum, CI must run formatting/linting, type-aware static checks,
  unit tests, CLI smoke tests, and `python -m compileall`.
- Mock external assemblers for unit tests. Keep a separately marked integration
  test for installed bioinformatics tools and a small end-to-end fixture.
- For algorithm changes, save machine-readable before/after metrics and the
  exact command. If expected scientific output changes, document why.

## 8. Documentation and repository hygiene

- Keep `README.md`, `Pipeline.md`, CLI help, and defaults synchronized.
- Document Python packages and external tools separately, including supported
  versions and how PHap locates them.
- Examples must use portable placeholder paths, distinct Hi-C mate files, and
  commands that work from a clean checkout.
- Do not commit generated assemblies, BAM/PAF intermediates, logs, caches,
  `.DS_Store`, `__pycache__`, or credentials. Add generated patterns to
  `.gitignore`.
- Before handoff, report files changed, checks run, checks not run, and any
  expected scientific-output difference.
