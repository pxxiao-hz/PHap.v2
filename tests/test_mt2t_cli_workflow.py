from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class Mt2tCliWorkflowTests(unittest.TestCase):
    def run_mt2t(self, fasta: Path, output: Path, min_length: int) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                "-m",
                "scripts.phap_mT2T",
                "--p_ctg",
                str(fasta),
                "--output-directory",
                str(output),
                "--min_contig_length",
                str(min_length),
                "--min_alignment_length",
                "1",
                "--threads",
                "1",
                "--process",
                "1",
                "--cpu-budget",
                "1",
            ],
            cwd=ROOT,
            env={**os.environ, "PYTHONPATH": str(ROOT)},
            text=True,
            capture_output=True,
            check=False,
        )

    def test_single_eligible_contig_bypasses_external_tools(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "input.fa"
            fasta.write_text(">only source header\nACGTACGT\n", encoding="utf-8")
            output = root / "out"
            completed = self.run_mt2t(fasta, output, min_length=1)
            assembled = (
                output / "04.remove.redundancy" / "mT2T.fa"
            ).read_text(encoding="utf-8")
            merge_paf = (output / "03.alignment" / "merge.paf").read_text(
                encoding="utf-8"
            )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(merge_paf, "")
        self.assertEqual(assembled, ">only source header\nACGTACGT\n")

    def test_zero_eligible_contigs_are_preserved_and_rerun_is_not_stale(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "input.fa"
            output = root / "out"
            fasta.write_text(">short\nAAAA\n", encoding="utf-8")
            first = self.run_mt2t(fasta, output, min_length=10)
            first_manifest = (output / "mt2t_upstream_manifest.tsv").read_text(
                encoding="utf-8"
            )
            fasta.write_text(">short\nCCCC\n", encoding="utf-8")
            second = self.run_mt2t(fasta, output, min_length=10)
            second_manifest = (output / "mt2t_upstream_manifest.tsv").read_text(
                encoding="utf-8"
            )
            assembled = (
                output / "04.remove.redundancy" / "mT2T.fa"
            ).read_text(encoding="utf-8")
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertNotEqual(first_manifest, second_manifest)
        self.assertEqual(assembled, ">short\nCCCC\n")

    def test_unsafe_identifier_is_never_used_as_a_path_or_shell_command(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "input.fa"
            fasta.write_text(">../unsafe;touch_BAD description\nACGT\n", encoding="utf-8")
            output = root / "out"
            completed = self.run_mt2t(fasta, output, min_length=1)
            id_map = (output / "mt2t_split_id_map.tsv").read_text(encoding="utf-8")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("contig_000001\t../unsafe;touch_BAD", id_map)

    def test_cpu_budget_is_validated_before_external_work(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "input.fa"
            fasta.write_text(">a\nAAAA\n>b\nCCCC\n", encoding="utf-8")
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "scripts.phap_mT2T",
                    "--p_ctg",
                    str(fasta),
                    "--output-directory",
                    str(root / "out"),
                    "--min_contig_length",
                    "1",
                    "--threads",
                    "2",
                    "--process",
                    "2",
                    "--cpu-budget",
                    "2",
                ],
                cwd=ROOT,
                env={**os.environ, "PYTHONPATH": str(ROOT)},
                text=True,
                capture_output=True,
                check=False,
            )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("exceeds --cpu-budget", completed.stderr)


if __name__ == "__main__":
    unittest.main()
