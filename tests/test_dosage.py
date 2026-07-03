from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from phap_core.dosage import (
    DosageError,
    WindowDepth,
    classify_window,
    classify_windows,
    fit_dosage_model,
    parse_window_depths,
    summarize_unitigs,
)


ROOT = Path(__file__).resolve().parents[1]
ENTRY_POINT = ROOT / "PHap.py"
LEGACY_ENTRY_POINT = ROOT / "shell" / "dosage.analysis.contig.type.identified.py"


class DosageModelTests(unittest.TestCase):
    def test_auto_fit_recovers_integer_copy_spacing(self) -> None:
        depths = [19.0, 20.0, 21.0, 39.0, 40.0, 41.0, 59.0, 60.0, 61.0, 79.0, 80.0, 81.0]
        model = fit_dosage_model(depths, ploidy=4)
        self.assertAlmostEqual(model.haploid_depth, 20.0, delta=1.0)
        self.assertEqual(model.estimation_method, "robust_integer_spacing")
        self.assertEqual(model.estimation_status, "well_separated")

    def test_explicit_haploid_depth_overrides_auto_fit(self) -> None:
        model = fit_dosage_model(
            [30.0, 60.0, 90.0],
            ploidy=3,
            haploid_depth=25.0,
            relative_sigma=0.1,
        )
        self.assertEqual(model.haploid_depth, 25.0)
        self.assertEqual(model.relative_sigma, 0.1)
        self.assertEqual(model.estimation_status, "user_supplied")

    def test_unidentifiable_single_peak_is_reported(self) -> None:
        model = fit_dosage_model([40.0] * 20, ploidy=4)
        self.assertEqual(model.estimation_status, "weakly_identified")
        self.assertTrue(model.alternatives)
        self.assertEqual(model.alternatives[0].score, model.objective_score)

    def test_window_states_include_low_ambiguous_and_high_copy(self) -> None:
        model = fit_dosage_model(
            [20.0, 40.0, 60.0, 80.0],
            ploidy=4,
            haploid_depth=20.0,
            relative_sigma=0.1,
        )
        low = classify_window(WindowDepth("low", 0, 10, 5.0), model)
        ambiguous = classify_window(WindowDepth("amb", 0, 10, 30.0), model)
        high = classify_window(WindowDepth("high", 0, 10, 100.1), model)
        self.assertEqual(low.classification, "low_coverage")
        self.assertEqual(ambiguous.classification, "ambiguous")
        self.assertEqual(high.classification, "high_copy")
        self.assertIsNone(ambiguous.dosage)

    def test_mixed_windows_are_not_hidden_by_unitig_average(self) -> None:
        model = fit_dosage_model(
            [20.0, 40.0],
            ploidy=4,
            haploid_depth=20.0,
            relative_sigma=0.1,
        )
        calls = classify_windows(
            [
                WindowDepth("mixed_utg", 10, 20, 40.0),
                WindowDepth("mixed_utg", 0, 10, 20.0),
            ],
            model,
        )
        summary = summarize_unitigs(calls)[0]
        self.assertAlmostEqual(summary.average_depth, 30.0)
        self.assertEqual(summary.status, "mixed")
        self.assertEqual(summary.contig_type, "ambiguous")
        self.assertTrue(summary.mixed_dosage)
        self.assertEqual(
            dict(summary.class_counts),
            {"dosage_1": 1, "dosage_2": 1},
        )

    def test_invalid_depth_is_rejected_with_line_number(self) -> None:
        with self.assertRaisesRegex(DosageError, "line 1"):
            parse_window_depths(
                ["utg1\t0\t100\t-1\n"],
                contig_column=0,
                start_column=1,
                end_column=2,
                depth_column=3,
            )


class DosageCliTests(unittest.TestCase):
    def _write_pandepth_fixture(self, directory: Path) -> Path:
        input_path = directory / "windows.tsv"
        rows = [
            ("utg_hap", 0, 10_000, 20.0),
            ("utg_hap", 10_000, 20_000, 20.0),
            ("utg_di", 0, 10_000, 40.0),
            ("utg_di", 10_000, 20_000, 40.0),
            ("utg_tri", 0, 10_000, 60.0),
            ("utg_tri", 10_000, 20_000, 60.0),
            ("utg_tetra", 0, 10_000, 80.0),
            ("utg_tetra", 10_000, 20_000, 80.0),
            ("utg_mixed", 0, 10_000, 20.0),
            ("utg_mixed", 10_000, 20_000, 40.0),
        ]
        with input_path.open("w", encoding="utf-8", newline="\n") as handle:
            for contig, start, end, depth in rows:
                handle.write(
                    f"{contig}\t{start}\t{end}\t10000\t.\t.\t.\t{depth}\n"
                )
        return input_path

    def test_installable_cli_writes_auditable_stable_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            input_path = self._write_pandepth_fixture(directory)
            command = [
                sys.executable,
                str(ENTRY_POINT),
                "dosage",
                "--input-file",
                str(input_path),
                "--pandepth",
                "--ploidy",
                "4",
                "--output-dir",
                str(directory),
            ]
            first = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
            self.assertEqual(first.returncode, 0, first.stderr)
            unitig_path = directory / "contig_depth.txt"
            window_path = directory / "dosage_windows.tsv"
            model_path = directory / "dosage_model.json"
            first_outputs = tuple(
                path.read_bytes() for path in (unitig_path, window_path, model_path)
            )

            second = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
            self.assertEqual(second.returncode, 0, second.stderr)
            second_outputs = tuple(
                path.read_bytes() for path in (unitig_path, window_path, model_path)
            )
            self.assertEqual(first_outputs, second_outputs)

            lines = unitig_path.read_text(encoding="utf-8").splitlines()
            self.assertEqual(
                lines[0].split("\t")[:3],
                ["contig_ID", "average_depth", "contig_type"],
            )
            mixed_fields = next(
                line.split("\t") for line in lines if line.startswith("utg_mixed\t")
            )
            self.assertEqual(mixed_fields[2], "ambiguous")
            self.assertEqual(mixed_fields[4], "mixed")
            model_payload = json.loads(model_path.read_text(encoding="utf-8"))
            self.assertEqual(model_payload["ploidy"], 4)
            self.assertIn("haploid_depth", model_payload)

    def test_legacy_script_path_remains_usable(self) -> None:
        result = subprocess.run(
            [sys.executable, str(LEGACY_ENTRY_POINT), "--help"],
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--ploidy", result.stdout)

    def test_cli_rejects_weakly_identified_auto_model(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            input_path = directory / "single_peak.tsv"
            input_path.write_text(
                "".join(
                    f"utg{index}\t0\t10000\t10000\t.\t.\t.\t40\n"
                    for index in range(4)
                ),
                encoding="utf-8",
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(ENTRY_POINT),
                    "dosage",
                    "--input-file",
                    str(input_path),
                    "--pandepth",
                    "--ploidy",
                    "4",
                    "--output-dir",
                    str(directory),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 1)
            self.assertIn("weakly identified", result.stderr)
            self.assertIn("--haploid-depth", result.stderr)
            self.assertFalse((directory / "contig_depth.txt").exists())


if __name__ == "__main__":
    unittest.main()
