from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class LocusWorkflowTests(unittest.TestCase):
    def test_complete_paf_to_locus_routing_and_union_bin_support(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            fasta = root / "unitigs.fa"
            fasta.write_text(
                f">u1 source description\n{'A' * 1000}\n"
                f">u2\n{'C' * 1000}\n",
                encoding="utf-8",
            )
            dosage = root / "dosage.tsv"
            dosage.write_text(
                "contig_ID\taverage_depth\tcontig_type\tdosage\tstatus\n"
                "u1\t20\thaplotig\t1\tassigned\n"
                "u2\t20\thaplotig\t1\tassigned\n",
                encoding="utf-8",
            )
            paf = root / "complete.paf"
            paf.write_text(
                "u1\t1000\t0\t600\t+\tchr1\t1000\t0\t600"
                "\t570\t600\t60\tcg:Z:600M\ttp:A:P\n"
                "u1\t1000\t400\t1000\t+\tchr1\t1000\t400\t1000"
                "\t570\t600\t60\ttp:A:P\tcg:Z:600M\n"
                "u2\t1000\t0\t900\t+\tchr1\t1000\t50\t950"
                "\t855\t900\t60\ttp:A:P\n",
                encoding="utf-8",
            )
            filtered = root / "filtered.paf"
            alignment_audit = root / "alignment.tsv"
            candidates = root / "candidates.tsv"
            decisions = root / "decisions.tsv"
            routing = root / "routing.tsv"
            manifest = root / "manifest.tsv"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "utils.paf_locus_evidence",
                    "--paf",
                    str(paf),
                    "--unitig-fasta",
                    str(fasta),
                    "--contig-type",
                    str(dosage),
                    "--target-count",
                    "1",
                    "--min-query-length",
                    "1",
                    "--min-alignment-block-length",
                    "1",
                    "--max-query-gap",
                    "1000",
                    "--max-target-gap",
                    "1000",
                    "--min-locus-identity",
                    "0.8",
                    "--min-locus-query-coverage",
                    "0.5",
                    "--min-locus-score-margin",
                    "0.05",
                    "--filtered-paf",
                    str(filtered),
                    "--alignment-audit",
                    str(alignment_audit),
                    "--candidate-audit",
                    str(candidates),
                    "--decision-audit",
                    str(decisions),
                    "--routing-audit",
                    str(routing),
                    "--manifest",
                    str(manifest),
                ],
                cwd=root,
                env={**os.environ, "PYTHONPATH": str(ROOT)},
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            top = root / "top.tsv"
            allelic = root / "allelic.tsv"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "utils.allelic_table_generate",
                    "--paf_file",
                    str(filtered),
                    "--min_align_length",
                    "1",
                    "--min_unitig_length",
                    "1",
                    "--bin_size",
                    "1000",
                    "--top_n",
                    "2",
                    "--chr_num",
                    "1",
                    "--contig_type",
                    str(dosage),
                    "--out_top_contigs_per_bin",
                    str(top),
                    "--out_allelic_table",
                    str(allelic),
                ],
                cwd=root,
                env={**os.environ, "PYTHONPATH": str(ROOT)},
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            sequence_directory = root / "sequences"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "utils.extract_locus_sequences",
                    "--p-utg",
                    str(fasta),
                    "--locus-decisions",
                    str(decisions),
                    "--output-directory",
                    str(sequence_directory),
                ],
                cwd=root,
                env={**os.environ, "PYTHONPATH": str(ROOT)},
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            top_text = top.read_text(encoding="utf-8")
            assigned = (sequence_directory / "chr1.putg.fa").read_text(
                encoding="utf-8"
            )
            sequence_routing = (
                sequence_directory / "locus_sequence_routing.tsv"
            ).read_text(encoding="utf-8")
            manifest_text = manifest.read_text(encoding="utf-8")
            alignment_audit_text = alignment_audit.read_text(encoding="utf-8")
            candidate_text = candidates.read_text(encoding="utf-8")

        self.assertIn("chr1\t0\t1000\tu1\t1000", top_text)
        self.assertIn("chr1\t0\t1000\tu2\t900", top_text)
        self.assertIn(">u1 source description", assigned)
        self.assertIn(">u2", assigned)
        self.assertIn("u1\tchr1.putg.fa", sequence_routing)
        self.assertIn("u2\tchr1.putg.fa", sequence_routing)
        self.assertIn("primary_policy\ttp:A:P_or_I", manifest_text)
        self.assertIn("allowed_loci\tchr1", manifest_text)
        self.assertIn("final_status\tfinal_reason", alignment_audit_text)
        self.assertIn("\twritten\tselected_assigned_locus", alignment_audit_text)
        self.assertIn("selected_for_output\tselection_reason", candidate_text)


if __name__ == "__main__":
    unittest.main()
