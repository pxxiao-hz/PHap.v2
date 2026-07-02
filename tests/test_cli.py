from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENTRY_POINT = ROOT / "PHap.py"


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(ENTRY_POINT), *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


class CliSmokeTests(unittest.TestCase):
    def test_top_level_help(self) -> None:
        result = run_cli("--help")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Usage: phap", result.stdout)
        self.assertIn("phase_reads", result.stdout)

    def test_version(self) -> None:
        result = run_cli("--version")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Version: 1.1.0", result.stdout)

    def test_subcommand_help_does_not_require_scientific_imports(self) -> None:
        for command in ("mt2t", "cluster", "phase_reads"):
            with self.subTest(command=command):
                result = run_cli(command, "--help")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("usage:", result.stdout.lower())

    def test_unknown_command_is_an_error(self) -> None:
        result = run_cli("unknown")
        self.assertEqual(result.returncode, 2)
        self.assertIn("unknown command", result.stderr)

    def test_preflight_reports_missing_input_without_traceback(self) -> None:
        result = run_cli("mt2t", "--p_ctg", "missing.fa")
        self.assertEqual(result.returncode, 1)
        self.assertIn("missing input file", result.stderr)
        self.assertNotIn("Traceback", result.stderr)


if __name__ == "__main__":
    unittest.main()
