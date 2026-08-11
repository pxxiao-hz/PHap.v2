from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
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
POLICY_EVIDENCE = ROOT / "tests" / "fixtures" / "dosage_policy_before_after.tsv"


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

    def test_window_states_include_low_depth_dosage_one_ambiguous_and_high_copy(
        self,
    ) -> None:
        model = fit_dosage_model(
            [20.0, 40.0, 60.0, 80.0],
            ploidy=4,
            haploid_depth=20.0,
            relative_sigma=0.1,
        )
        low = classify_window(WindowDepth("low", 0, 10, 5.0), model)
        ambiguous = classify_window(WindowDepth("amb", 0, 10, 30.0), model)
        high = classify_window(WindowDepth("high", 0, 10, 100.1), model)
        self.assertEqual(low.classification, "dosage_1")
        self.assertEqual(low.dosage, 1)
        self.assertEqual(ambiguous.classification, "ambiguous")
        self.assertEqual(high.classification, "high_copy")
        self.assertIsNone(ambiguous.dosage)

    def test_unitig_call_uses_average_depth_not_window_class_fraction(self) -> None:
        model = fit_dosage_model(
            [20.0, 40.0],
            ploidy=4,
            haploid_depth=20.0,
            relative_sigma=0.1,
        )
        calls = classify_windows(
            [
                WindowDepth("mixed_utg", 0, 10, 20.0),
                WindowDepth("mixed_utg", 10, 20, 20.0),
                WindowDepth("mixed_utg", 20, 30, 20.0),
                WindowDepth("mixed_utg", 30, 40, 30.0),
                WindowDepth("mixed_utg", 40, 50, 30.0),
            ],
            model,
        )
        summary = summarize_unitigs(calls, model)[0]
        self.assertAlmostEqual(summary.average_depth, 24.0)
        self.assertEqual(summary.status, "assigned")
        self.assertEqual(summary.contig_type, "haplotig")
        self.assertEqual(summary.dosage, 1)
        self.assertEqual(summary.dominant_class, "dosage_1")
        self.assertEqual(summary.dominant_fraction, 0.6)
        self.assertFalse(summary.mixed_dosage)
        self.assertEqual(
            dict(summary.class_counts),
            {"ambiguous": 2, "dosage_1": 3},
        )

    def test_ambiguous_average_depth_remains_ambiguous(self) -> None:
        model = fit_dosage_model(
            [20.0, 40.0],
            ploidy=4,
            haploid_depth=20.0,
            relative_sigma=0.1,
        )
        calls = classify_windows(
            [
                WindowDepth("balanced", 0, 10, 20.0),
                WindowDepth("balanced", 10, 20, 40.0),
            ],
            model,
        )
        summary = summarize_unitigs(calls, model)[0]
        self.assertEqual(summary.status, "ambiguous")
        self.assertEqual(summary.contig_type, "ambiguous")
        self.assertIsNone(summary.dosage)
        self.assertTrue(summary.mixed_dosage)

    def test_average_depth_can_override_dominant_window_class(self) -> None:
        model = fit_dosage_model(
            [20.0, 40.0],
            ploidy=4,
            haploid_depth=20.0,
            relative_sigma=0.1,
        )
        calls = classify_windows(
            [
                WindowDepth("average_only", 0, 10, 20.0),
                WindowDepth("average_only", 10, 20, 20.0),
                WindowDepth("average_only", 20, 30, 20.0),
                WindowDepth("average_only", 30, 40, 60.0),
                WindowDepth("average_only", 40, 50, 60.0),
            ],
            model,
        )
        summary = summarize_unitigs(calls, model)[0]
        self.assertAlmostEqual(summary.average_depth, 36.0)
        self.assertEqual(summary.dominant_class, "dosage_1")
        self.assertEqual(summary.dominant_fraction, 0.6)
        self.assertEqual(summary.dosage, 2)
        self.assertEqual(summary.contig_type, "diplotig")
        self.assertEqual(summary.status, "assigned")
        self.assertTrue(summary.mixed_dosage)

    def test_machine_readable_policy_evidence_matches_current_calls(self) -> None:
        with POLICY_EVIDENCE.open(encoding="utf-8", newline="") as handle:
            rows = tuple(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(rows)
        for row in rows:
            with self.subTest(case_id=row["case_id"]):
                depths = [float(value) for value in row["depths"].split(",")]
                model = fit_dosage_model(
                    depths,
                    ploidy=4,
                    haploid_depth=float(row["haploid_depth"]),
                    relative_sigma=0.1,
                )
                calls = classify_windows(
                    [
                        WindowDepth(row["case_id"], index, index + 1, depth)
                        for index, depth in enumerate(depths)
                    ],
                    model,
                )
                summary = summarize_unitigs(
                    calls,
                    model,
                    min_confidence=float(row["min_confidence"]),
                )[0]
                observed_counts = ",".join(
                    f"{name}:{count}"
                    for name, count in sorted(
                        Counter(call.classification for call in calls).items()
                    )
                )
                self.assertEqual(observed_counts, row["class_counts"])
                self.assertAlmostEqual(
                    summary.average_depth,
                    float(row["average_depth"]),
                )
                self.assertEqual(summary.status, row["new_status"])
                self.assertEqual(summary.contig_type, row["new_contig_type"])
                observed_dosage = "." if summary.dosage is None else str(summary.dosage)
                self.assertEqual(observed_dosage, row["new_dosage"])

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
            ("utg_dominant", 0, 10_000, 20.0),
            ("utg_dominant", 10_000, 20_000, 20.0),
            ("utg_dominant", 20_000, 30_000, 20.0),
            ("utg_dominant", 30_000, 40_000, 30.0),
            ("utg_dominant", 40_000, 50_000, 30.0),
            ("utg_balanced", 0, 10_000, 20.0),
            ("utg_balanced", 10_000, 20_000, 40.0),
            ("utg_low", 0, 10_000, 5.0),
            ("utg_low", 10_000, 20_000, 5.0),
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
            rows_by_id = {
                fields[0]: fields
                for fields in (line.split("\t") for line in lines[1:])
            }
            dominant_fields = rows_by_id["utg_dominant"]
            self.assertEqual(dominant_fields[2:7], [
                "haplotig",
                "1",
                "assigned",
                "dosage_1",
                "0.600000",
            ])
            self.assertEqual(dominant_fields[7], "false")
            balanced_fields = rows_by_id["utg_balanced"]
            self.assertEqual(balanced_fields[2], "ambiguous")
            self.assertEqual(balanced_fields[4], "ambiguous")
            low_fields = rows_by_id["utg_low"]
            self.assertEqual(low_fields[2:5], ["haplotig", "1", "assigned"])

            low_window_fields = next(
                line.split("\t")
                for line in window_path.read_text(encoding="utf-8").splitlines()
                if line.startswith("utg_low\t")
            )
            self.assertEqual(low_window_fields[5:7], ["dosage_1", "1"])
            model_payload = json.loads(model_path.read_text(encoding="utf-8"))
            self.assertEqual(model_payload["ploidy"], 4)
            self.assertIn("haploid_depth", model_payload)
            self.assertEqual(
                model_payload["classification"]["low_depth_policy"],
                "dosage_1",
            )
            self.assertEqual(
                model_payload["classification"]["unitig_summary_policy"],
                "average_depth",
            )
            self.assertNotIn(
                "min_unitig_support",
                model_payload["classification"],
            )

    def test_removed_unitig_support_option_is_rejected(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(ENTRY_POINT),
                "dosage",
                "--input-file",
                "unused.tsv",
                "--ploidy",
                "4",
                "--min-unitig-support",
                "0.5",
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("unrecognized arguments: --min-unitig-support", result.stderr)

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
