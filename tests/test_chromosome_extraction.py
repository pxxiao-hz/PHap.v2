import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXTRACT_SCRIPT = ROOT / "utils" / "extract_chr_from_putg.py"


def paf_line(
    query,
    query_length,
    query_start,
    query_end,
    target,
    target_length,
    target_start,
    target_end,
    tags,
    strand="+",
):
    block_length = query_end - query_start
    fields = [
        query,
        query_length,
        query_start,
        query_end,
        strand,
        target,
        target_length,
        target_start,
        target_end,
        block_length,
        block_length,
        60,
        *tags,
    ]
    return "\t".join(map(str, fields))


class ChromosomeExtractionTests(unittest.TestCase):
    def run_extractor(self, work, fasta_text, paf_text, check=True, chain_qc_text=None):
        p_utg = work / "p_utg.fa"
        reference = work / "mT2T.fa"
        paf = work / "chain.paf"
        output = work / "02.chr_seq"
        p_utg.write_text(fasta_text)
        reference.write_text(">chr1\n" + "A" * 100 + "\n>chr2\n" + "C" * 80 + "\n")
        paf.write_text(paf_text)
        command = [
            sys.executable,
            str(EXTRACT_SCRIPT),
            "--p_utg",
            str(p_utg),
            "--mT2T",
            str(reference),
            "--paf",
            str(paf),
            "--wd",
            str(output),
            "--chr_num",
            "2",
        ]
        if chain_qc_text is not None:
            chain_qc = work / "chain.qc.tsv"
            chain_qc.write_text(chain_qc_text)
            command.extend(["--chain-qc", str(chain_qc)])
        result = subprocess.run(
            command,
            check=check,
            capture_output=True,
            text=True,
        )
        return output, result

    def test_assigns_selected_chain_and_handles_tags_by_name(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            output = work / "02.chr_seq"
            output.mkdir()
            (output / "chr99.putg.fa").write_text(">stale\nA\n")
            (output / "contig_match_ratios.xlsx").write_text("stale")

            records = [
                # Overlapping query intervals must have 12 bp union, not 14 bp sum.
                paf_line(
                    "q1", 12, 0, 8, "chr1", 100, 0, 8,
                    ["tp:A:P", "cg:Z:8M", "pv:Z:dense"],
                ),
                paf_line(
                    "q1", 12, 6, 12, "chr1", 100, 6, 12,
                    ["pv:Z:dense", "tp:A:P", "cg:Z:6M"],
                ),
                paf_line(
                    "q2", 8, 0, 8, "chr1", 100, 20, 28,
                    ["cg:Z:8M", "pv:Z:dense", "tp:A:S"],
                ),
                paf_line(
                    "q3", 10, 0, 10, "chr1", 100, 30, 40,
                    ["tp:A:P", "pv:Z:dense"],
                ),
                paf_line(
                    "q3", 10, 0, 10, "chr2", 80, 30, 40,
                    ["tp:A:P", "pv:Z:dense"],
                ),
            ]
            fasta = ">q1 retained description\nAAAAAAAAAAAA\n>q2\nCCCCCCCC\n>q3\nGGGGGGGGGG\n>q4\nTTTT\n"
            output, result = self.run_extractor(work, fasta, "\n".join(records) + "\n")

            self.assertIn("1 assigned, 3 unassigned", result.stdout)
            self.assertEqual(
                (output / "chr1.putg.fa").read_text(),
                ">q1 retained description\nAAAAAAAAAAAA\n",
            )
            self.assertEqual(
                set(line[1:].split()[0] for line in (output / "un_chr.fa").read_text().splitlines() if line.startswith(">")),
                {"q2", "q3", "q4"},
            )
            self.assertFalse((output / "chr99.putg.fa").exists())
            self.assertFalse((output / "contig_match_ratios.xlsx").exists())

            with (output / "chromosome_assignments.tsv").open() as handle:
                rows = {row["unitig"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["q1"]["query_aligned_bp"], "12")
            self.assertEqual(rows["q1"]["query_coverage"], "1.000000")
            self.assertEqual(rows["q1"]["placement_mode"], "dense")
            self.assertEqual(rows["q2"]["reason"], "secondary_only")
            self.assertEqual(rows["q3"]["reason"], "ambiguous_target")
            self.assertEqual(rows["q4"]["reason"], "no_accepted_collinear_chain")

            summary = json.loads((output / "chromosome_assignment.summary.json").read_text())
            self.assertEqual(summary["assigned"]["unitigs"], 1)
            self.assertEqual(summary["unassigned"]["unitigs"], 3)
            self.assertEqual(summary["chromosome_selection_source"], "reference_fasta")

    def test_accepts_mixed_strands_when_the_selected_target_is_unique(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            records = [
                paf_line(
                    "q1", 12, 0, 6, "chr1", 100, 0, 6,
                    ["tp:A:P", "pv:Z:segmented"],
                    strand="+",
                ),
                paf_line(
                    "q1", 12, 6, 12, "chr1", 100, 20, 26,
                    ["tp:A:P", "pv:Z:segmented"],
                    strand="-",
                ),
            ]
            chain_qc = (
                "query\tquery_length\tbest_target\tbest_strand\treference_margin\tstatus\treason\n"
                "q1\t12\tchr1\tmixed\t1.0\taccepted\t\n"
            )
            output, _ = self.run_extractor(
                work,
                ">q1\nAAAAAAAAAAAA\n",
                "\n".join(records) + "\n",
                chain_qc_text=chain_qc,
            )

            with (output / "chromosome_assignments.tsv").open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["status"], "assigned")
            self.assertEqual(row["target"], "chr1")
            self.assertEqual(row["strand"], "mixed")
            self.assertEqual(row["placement_mode"], "segmented")

    def test_input_mismatch_does_not_replace_existing_outputs(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            output = work / "02.chr_seq"
            output.mkdir()
            existing = output / "chr1.putg.fa"
            existing.write_text(">old\nAAAA\n")
            record = paf_line(
                "q1", 12, 0, 8, "chr1", 100, 0, 8,
                ["tp:A:P", "cg:Z:8M", "pv:Z:dense"],
            )
            _, result = self.run_extractor(work, ">q1\nAAAAAAAAAAA\n", record + "\n", check=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("does not match the p_utg FASTA", result.stderr)
            self.assertEqual(existing.read_text(), ">old\nAAAA\n")
            self.assertEqual(list(output.glob(".chr_seq.*")), [])

    def test_rescues_clear_best_target_but_retains_ambiguous_chain(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            record = paf_line(
                "q1", 12, 0, 12, "chr1", 100, 0, 12,
                ["tp:A:P", "pv:Z:dense"],
            )
            chain_qc = (
                "query\tquery_length\tbest_target\tbest_strand\treference_margin\tstatus\treason\n"
                "q1\t12\tchr1\t+\t0.90\taccepted\t\n"
                "q2\t8\tchr2\t+\t0.80\trejected\tquery_coverage_below_threshold\n"
                "q3\t6\tchr2\t+\t0.01\trejected\tchain_aligned_bp_below_threshold\n"
            )
            output, _ = self.run_extractor(
                work,
                ">q1\nAAAAAAAAAAAA\n>q2\nCCCCCCCC\n>q3\nGGGGGG\n",
                record + "\n",
                chain_qc_text=chain_qc,
            )
            with (output / "chromosome_assignments.tsv").open() as handle:
                rows = {row["unitig"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(
                rows["q2"]["reason"],
                "rescued_best_chromosome:query_coverage_below_threshold",
            )
            self.assertEqual(rows["q2"]["status"], "assigned")
            self.assertEqual(rows["q2"]["target"], "chr2")
            self.assertEqual(rows["q2"]["reference_margin"], "0.800000")
            self.assertEqual(rows["q3"]["status"], "unassigned")
            self.assertEqual(
                rows["q3"]["reason"],
                "chromosome_margin_below_threshold:chain_aligned_bp_below_threshold",
            )


if __name__ == "__main__":
    unittest.main()
