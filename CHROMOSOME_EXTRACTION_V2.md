# PHap v2 Chromosome Sequence Extraction

## Assignment Rule

`02.chr_seq` consumes `putg_vs_mT2T.collinear.paf` and
`collinear_chain.qc.tsv`. The preceding selection step records the best
reference target from aggregate local evidence. Its strand can be `mixed`
because inversions and direction changes within one unitig are valid.

Chromosome membership deliberately uses a less restrictive decision than the
allelic table. An accepted placement is assigned directly. A rejected placement
is rescued when it has a reference target and `(best_score - second_score) /
best_score` is at least `--min-chromosome-margin` (default 0.05). Query
coverage, identity, aligned length, and target-span coverage are retained
as QC metrics but do not block chromosome assignment. Unitigs without a target
or with insufficient best-target separation remain unassigned.

The mT2T FASTA defines the chromosome set independently of which targets have
PAF hits. The longest `--chr_num` reference records are selected. PAF optional
tags such as `tp`, `cg`, and `pv` are parsed by tag name, never by column
position.

Assigned and rescued unitigs are written to `<target>.putg.fa`. Ambiguous and
unaligned unitigs are written once to `un_chr.fa`. Sequence records are
streamed, so memory use does not scale with total p_utg FASTA size.

## Safety And Consistency Checks

- PAF query and target coordinates must be valid.
- Query lengths must agree across PAF, chain QC, and p_utg FASTA.
- PAF target lengths must agree with the mT2T FASTA.
- More than one selected target for a unitig is rejected as ambiguous; both
  strands on the same selected target are combined and reported as `mixed`.
- Low allelic-projection coverage does not discard a unitig whose chromosome
  target is unambiguous.
- Overlapping query intervals are unioned for QC coverage instead of counted
  more than once.
- Duplicate FASTA IDs and PAF IDs absent from p_utg are fatal errors.
- New outputs are built in a temporary directory. Existing results are only
  replaced after all validation succeeds.
- Stale `*.putg.fa` files from earlier runs and obsolete v1 matrix files are
  removed after successful generation.

## Outputs

| File | Meaning |
| --- | --- |
| `chrXX.putg.fa` | Unitigs assigned to one accepted chromosome chain |
| `un_chr.fa` | Rejected or unaligned unitigs retained for the rescue step |
| `chromosome_assignments.tsv` | One row per input unitig with assignment/rescue reason, placement mode, score margin, and chain metrics |
| `chromosome_assignment.summary.json` | Per-chromosome counts/bp, rejection reasons, inputs, and conservation totals |
| `putg_vs_mT2T.best.match.txt` | Two-column compatibility view of accepted assignments |

The extractor does not use the allelic table or GFA. Those inputs determine
high-confidence mutual-exclusion constraints for haplotype clustering;
chromosome membership only requires an adequately separated best target.
