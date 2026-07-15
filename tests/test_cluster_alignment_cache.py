from __future__ import annotations

import sys
import tempfile
import threading
import time
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Sequence, TextIO

from phap_core.stage_manifest import CacheValidation
from utils.cluster_alignment_workflow import build_or_reuse_locus_alignment


class FakeAlignmentRunner:
    def __init__(self, text: str = "alignment\n") -> None:
        self.text = text
        self.calls = 0
        self._lock = threading.Lock()

    def __call__(self, command: Sequence[str], *, stdout: TextIO) -> None:
        with self._lock:
            self.calls += 1
        stdout.write(self.text)
        return None


class ClusterAlignmentCacheTests(unittest.TestCase):
    def prepare(self, root: Path) -> tuple[Path, Path, Path, Path]:
        mt2t = root / "mt2t.fa"
        unitigs = root / "unitigs.fa"
        paf = root / "alignment.paf"
        manifest = root / "alignment.manifest.json"
        mt2t.write_text(">chr1\nAAAA\n", encoding="utf-8")
        unitigs.write_text(">u1\nCCCC\n", encoding="utf-8")
        return mt2t, unitigs, paf, manifest

    def run_stage(
        self,
        mt2t: Path,
        unitigs: Path,
        paf: Path,
        manifest: Path,
        runner: Any,
        *,
        threads: int = 2,
        tool_version: str = "2.30",
        phap_version: str = "test",
    ) -> CacheValidation:
        return build_or_reuse_locus_alignment(
            mt2t_fasta=str(mt2t),
            unitig_fasta=str(unitigs),
            output_paf=str(paf),
            manifest_path=str(manifest),
            threads=threads,
            minimap2_version=tool_version,
            minimap2_executable=sys.executable,
            phap_version=phap_version,
            runner=runner,
        )

    def test_identical_request_hits_without_running_alignment_twice(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            paths = self.prepare(Path(directory))
            runner = FakeAlignmentRunner()
            first = self.run_stage(*paths, runner)
            second = self.run_stage(*paths, runner)

        self.assertFalse(first.hit)
        self.assertTrue(second.hit)
        self.assertEqual(runner.calls, 1)

    def test_input_parameter_tool_and_phap_changes_invalidate_cache(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            mt2t, unitigs, paf, manifest = self.prepare(Path(directory))
            runner = FakeAlignmentRunner()
            self.run_stage(mt2t, unitigs, paf, manifest, runner)
            unitigs.write_text(">u1\nGGGG\n", encoding="utf-8")
            self.run_stage(mt2t, unitigs, paf, manifest, runner)
            self.run_stage(mt2t, unitigs, paf, manifest, runner, threads=3)
            self.run_stage(
                mt2t, unitigs, paf, manifest, runner, threads=3, tool_version="2.31"
            )
            self.run_stage(
                mt2t,
                unitigs,
                paf,
                manifest,
                runner,
                threads=3,
                tool_version="2.31",
                phap_version="next",
            )

        self.assertEqual(runner.calls, 5)

    def test_output_hash_corruption_and_orphan_output_force_rebuild(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            mt2t, unitigs, paf, manifest = self.prepare(Path(directory))
            runner = FakeAlignmentRunner("AAAA\n")
            self.run_stage(mt2t, unitigs, paf, manifest, runner)
            paf.write_text("BBBB\n", encoding="utf-8")
            self.run_stage(mt2t, unitigs, paf, manifest, runner)
            manifest.unlink()
            self.run_stage(mt2t, unitigs, paf, manifest, runner)

        self.assertEqual(runner.calls, 3)

    def test_failed_rebuild_preserves_previous_output_without_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            mt2t, unitigs, paf, manifest = self.prepare(Path(directory))
            good = FakeAlignmentRunner("complete\n")
            self.run_stage(mt2t, unitigs, paf, manifest, good)
            unitigs.write_text(">u1\nGGGG\n", encoding="utf-8")

            def failing_runner(command: Sequence[str], *, stdout: TextIO) -> None:
                stdout.write("partial\n")
                raise RuntimeError("alignment failed")

            with self.assertRaisesRegex(RuntimeError, "alignment failed"):
                self.run_stage(mt2t, unitigs, paf, manifest, failing_runner)
            observed = paf.read_text(encoding="utf-8")
            manifest_exists = manifest.exists()

        self.assertEqual(observed, "complete\n")
        self.assertFalse(manifest_exists)

    def test_empty_paf_is_a_valid_cached_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            paths = self.prepare(Path(directory))
            runner = FakeAlignmentRunner("")
            self.run_stage(*paths, runner)
            validation = self.run_stage(*paths, runner)

        self.assertTrue(validation.hit)
        self.assertEqual(runner.calls, 1)

    def test_input_change_during_alignment_prevents_manifest_commit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            mt2t, unitigs, paf, manifest = self.prepare(Path(directory))

            def mutating_runner(command: Sequence[str], *, stdout: TextIO) -> None:
                stdout.write("mixed\n")
                unitigs.write_text(">u1\nGGGG\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "inputs changed"):
                self.run_stage(mt2t, unitigs, paf, manifest, mutating_runner)
            manifest_exists = manifest.exists()

        self.assertFalse(manifest_exists)

    def test_concurrent_requests_run_alignment_once(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            paths = self.prepare(Path(directory))
            calls = 0
            calls_lock = threading.Lock()

            def slow_runner(command: Sequence[str], *, stdout: TextIO) -> None:
                nonlocal calls
                with calls_lock:
                    calls += 1
                time.sleep(0.05)
                stdout.write("alignment\n")

            with ThreadPoolExecutor(max_workers=2) as executor:
                futures = [
                    executor.submit(self.run_stage, *paths, slow_runner)
                    for _ in range(2)
                ]
                results = [future.result() for future in futures]

        self.assertEqual(calls, 1)
        self.assertEqual(sorted(result.hit for result in results), [False, True])

    def test_path_conflicts_are_rejected_before_any_input_is_modified(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            mt2t, unitigs, paf, manifest = self.prepare(Path(directory))
            original_mt2t = mt2t.read_bytes()
            original_unitigs = unitigs.read_bytes()
            runner = FakeAlignmentRunner()
            conflicts = (
                (mt2t, manifest),
                (paf, unitigs),
                (paf, paf),
            )
            for output, marker in conflicts:
                with self.subTest(output=output, marker=marker):
                    with self.assertRaisesRegex(ValueError, "paths must be distinct"):
                        self.run_stage(
                            mt2t,
                            unitigs,
                            output,
                            marker,
                            runner,
                        )

            self.assertEqual(mt2t.read_bytes(), original_mt2t)
            self.assertEqual(unitigs.read_bytes(), original_unitigs)
            self.assertEqual(runner.calls, 0)


if __name__ == "__main__":
    unittest.main()
