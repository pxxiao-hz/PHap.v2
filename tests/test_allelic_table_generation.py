from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class AllelicTableGenerationTests(unittest.TestCase):
    def run_generator(
        self,
        root: Path,
        paf_lines: list[str],
        dosage_rows: list[tuple[str, str, str]],
        *,
        ploidy: int = 4,
        chr_num: int = 1,
        legacy_compatible: bool = False,
        min_bin_support_bases: int = 1,
        min_bin_coverage: float = 0.0,
    ) -> tuple[subprocess.CompletedProcess[str], dict[str, Path]]:
        root.mkdir(parents=True, exist_ok=True)
        paf = root / "input.paf"
        paf.write_text("".join(paf_lines), encoding="utf-8")
        dosage = root / "dosage.tsv"
        dosage.write_text(
            "contig_ID\taverage_depth\tcontig_type\tdosage\tstatus\n"
            + "".join(
                f"{unitig}\t20\t{state}\t{copy_count}\tassigned\n"
                for unitig, state, copy_count in dosage_rows
            ),
            encoding="utf-8",
        )
        paths = {
            "top": root / "top.tsv",
            "allelic": root / "allelic.tsv",
            "candidates": root / "candidates.tsv",
            "bins": root / "bins.tsv",
        }
        command = [
                sys.executable,
                "-m",
                "utils.allelic_table_generate",
                "--paf_file",
                str(paf),
                "--contig_type",
                str(dosage),
                "--min_align_length",
                "1",
                "--min_unitig_length",
                "1",
                "--bin_size",
                "100",
                "--chr_num",
                str(chr_num),
                "--out_top_contigs_per_bin",
                str(paths["top"]),
                "--out_allelic_table",
                str(paths["allelic"]),
        ]
        if legacy_compatible:
            command.extend(["--top_n", str(ploidy)])
        else:
            command.extend(
                [
                    "--ploidy",
                    str(ploidy),
                    "--min-bin-support-bases",
                    str(min_bin_support_bases),
                    "--min-bin-coverage",
                    str(min_bin_coverage),
                    "--selection-audit",
                    str(paths["candidates"]),
                    "--bin-audit",
                    str(paths["bins"]),
                ]
            )
        completed = subprocess.run(
            command,
            cwd=root,
            env={**os.environ, "PYTHONPATH": str(ROOT)},
            text=True,
            capture_output=True,
            check=False,
        )
        return completed, paths

    def test_legacy_cli_reproduces_the_scientific_capacity_regression(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            completed, paths = self.run_generator(
                Path(directory),
                [
                    _paf("d3", "chr1", 0, 100),
                    _paf("d2", "chr1", 0, 99),
                    _paf("d1", "chr1", 0, 98),
                    _paf("x", "chr1", 0, 97),
                ],
                [
                    ("d3", "triplotig", "3"),
                    ("d2", "diplotig", "2"),
                    ("d1", "haplotig", "1"),
                    ("x", "haplotig", "1"),
                ],
                legacy_compatible=True,
            )
            allelic = paths["allelic"].read_text(encoding="utf-8")

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(allelic, "chr1\t0\t100\td3\td1\n")

    def test_capacity_safe_combination_is_selected_and_fully_audited(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paf_lines = [
                _paf("d3", "chr1", 0, 100),
                _paf("d2", "chr1", 0, 99),
                _paf("d1", "chr1", 0, 98),
                _paf("x", "chr1", 0, 97),
            ]
            completed, paths = self.run_generator(
                root,
                paf_lines,
                [
                    ("d3", "triplotig", "3"),
                    ("d2", "diplotig", "2"),
                    ("d1", "haplotig", "1"),
                    ("x", "haplotig", "1"),
                ],
            )
            top = paths["top"].read_text(encoding="utf-8")
            allelic = paths["allelic"].read_text(encoding="utf-8")
            candidate_rows = paths["candidates"].read_text(encoding="utf-8").splitlines()
            bin_audit = paths["bins"].read_text(encoding="utf-8")

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("chr1\t0\t100\td3\t100", top)
        self.assertIn("chr1\t0\t100\td1\t98", top)
        self.assertNotIn("\td2\t", top)
        self.assertEqual(allelic, "chr1\t0\t100\td3\td1\n")
        self.assertEqual(len(candidate_rows), 5)
        self.assertIn("copy_weighted_target_union_bases\t398\t1", bin_audit)

    def test_equal_optima_are_ambiguous_and_not_written_to_allelic_table(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            unitigs = ["collapsed", "h1", "h2", "h3", "h4"]
            completed, paths = self.run_generator(
                root,
                [_paf(unitig, "chr1", 0, 100) for unitig in unitigs],
                [("collapsed", "tetraplotig", "4")]
                + [(unitig, "haplotig", "1") for unitig in unitigs[1:]],
            )
            allelic = paths["allelic"].read_text(encoding="utf-8")
            candidate_rows = paths["candidates"].read_text(
                encoding="utf-8"
            ).splitlines()[1:]
            bin_audit = paths["bins"].read_text(encoding="utf-8")

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(allelic, "")
        self.assertTrue(all(row.split("\t")[9] == "ambiguous" for row in candidate_rows))
        self.assertIn("\tambiguous\tmultiple_equal_optimum_subsets\t", bin_audit)

    def test_multi_target_unitig_is_rejected_instead_of_order_tiebroken(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            completed, _ = self.run_generator(
                root,
                [_paf("u1", "chr2", 0, 100), _paf("u1", "chr1", 0, 100)],
                [("u1", "haplotig", "1")],
                chr_num=2,
            )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("assigned to multiple targets", completed.stderr)

    def test_bin_support_threshold_is_applied_and_audited(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            completed, paths = self.run_generator(
                Path(directory),
                [_paf("sliver", "chr1", 0, 1)],
                [("sliver", "haplotig", "1")],
                min_bin_support_bases=2,
            )
            allelic = paths["allelic"].read_text(encoding="utf-8")
            candidate_audit = paths["candidates"].read_text(encoding="utf-8")
            bin_audit = paths["bins"].read_text(encoding="utf-8")

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(allelic, "")
        self.assertIn("\tunassigned\tinsufficient_support_bases\t", candidate_audit)
        self.assertIn("\t2\t0.000000\tcopy_weighted_target_union_bases\t", bin_audit)

    def test_paf_input_order_does_not_change_outputs(self) -> None:
        paf_lines = [
            _paf("d2", "chr1", 0, 90),
            _paf("h1", "chr1", 0, 100),
            _paf("h2", "chr1", 0, 80),
        ]
        dosage = [
            ("d2", "diplotig", "2"),
            ("h1", "haplotig", "1"),
            ("h2", "haplotig", "1"),
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first_completed, first_paths = self.run_generator(
                root / "first", paf_lines, dosage
            )
            second_completed, second_paths = self.run_generator(
                root / "second", list(reversed(paf_lines)), dosage
            )
            first = {
                name: path.read_bytes() for name, path in first_paths.items()
            }
            second = {
                name: path.read_bytes() for name, path in second_paths.items()
            }

        self.assertEqual(first_completed.returncode, 0, first_completed.stderr)
        self.assertEqual(second_completed.returncode, 0, second_completed.stderr)
        self.assertEqual(first, second)


def _paf(unitig: str, target: str, start: int, end: int) -> str:
    length = end - start
    return (
        f"{unitig}\t100\t0\t{length}\t+\t{target}\t100\t{start}\t{end}\t"
        f"{length}\t{length}\t60\ttp:A:P\n"
    )


if __name__ == "__main__":
    unittest.main()
