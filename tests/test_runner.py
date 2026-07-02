from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path

from phap_core.runner import PreflightError, require_input_files, run_shell_commands_parallel


class RunnerTests(unittest.TestCase):
    def test_parallel_worker_failure_is_propagated(self) -> None:
        with self.assertRaises(subprocess.CalledProcessError):
            run_shell_commands_parallel(["exit 7"], max_workers=1)

    def test_missing_input_is_reported(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            missing = Path(directory) / "missing.fa"
            with self.assertRaisesRegex(PreflightError, "missing.fa"):
                require_input_files([missing])


if __name__ == "__main__":
    unittest.main()
