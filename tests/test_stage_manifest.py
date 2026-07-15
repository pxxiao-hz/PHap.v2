from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from phap_core.atomic_io import atomic_text_writer
from phap_core.stage_manifest import (
    StageInput,
    StageSpec,
    capture_stage,
    invalidate_stage_manifest,
    validate_stage_cache,
    write_stage_manifest,
)


class StageManifestTests(unittest.TestCase):
    def make_spec(
        self,
        source: Path,
        *,
        parameter: str = "asm5",
        tool_version: str = "2.30",
        phap_version: str = "test-version",
    ) -> StageSpec:
        return StageSpec(
            stage="test_alignment",
            inputs=(StageInput("query", source),),
            parameters=(("preset", parameter),),
            tool_versions=(("minimap2", tool_version),),
            phap_version=phap_version,
        )

    def test_matching_manifest_and_output_hash_are_a_cache_hit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.fa"
            output = root / "output.paf"
            manifest = root / "manifest.json"
            source.write_text(">u1\nAAAA\n", encoding="utf-8")
            output.write_text("alignment\n", encoding="utf-8")
            spec = self.make_spec(source)
            snapshot = capture_stage(spec)
            write_stage_manifest(snapshot, manifest, {"paf": output})
            validation = validate_stage_cache(snapshot, manifest, {"paf": output})

        self.assertTrue(validation.hit)
        self.assertEqual(validation.reason, "cache_valid")

    def test_input_content_change_invalidates_same_path(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.fa"
            output = root / "output.paf"
            manifest = root / "manifest.json"
            source.write_text(">u1\nAAAA\n", encoding="utf-8")
            output.write_text("alignment\n", encoding="utf-8")
            spec = self.make_spec(source)
            snapshot = capture_stage(spec)
            write_stage_manifest(snapshot, manifest, {"paf": output})
            source.write_text(">u1\nCCCC\n", encoding="utf-8")
            validation = validate_stage_cache(
                capture_stage(spec), manifest, {"paf": output}
            )

        self.assertFalse(validation.hit)
        self.assertEqual(validation.reason, "stage_signature_mismatch")

    def test_parameter_tool_and_phap_versions_are_part_of_signature(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.fa"
            output = root / "output.paf"
            manifest = root / "manifest.json"
            source.write_text(">u1\nAAAA\n", encoding="utf-8")
            output.write_text("alignment\n", encoding="utf-8")
            write_stage_manifest(
                capture_stage(self.make_spec(source)), manifest, {"paf": output}
            )
            validations = (
                validate_stage_cache(
                    capture_stage(self.make_spec(source, parameter="map-hifi")),
                    manifest,
                    {"paf": output},
                ),
                validate_stage_cache(
                    capture_stage(self.make_spec(source, tool_version="2.31")),
                    manifest,
                    {"paf": output},
                ),
                validate_stage_cache(
                    capture_stage(self.make_spec(source, phap_version="next")),
                    manifest,
                    {"paf": output},
                ),
            )

        self.assertTrue(all(not validation.hit for validation in validations))
        self.assertTrue(
            all(
                validation.reason == "stage_signature_mismatch"
                for validation in validations
            )
        )

    def test_missing_or_corrupt_output_is_never_a_hit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.fa"
            output = root / "output.paf"
            manifest = root / "manifest.json"
            source.write_text(">u1\nAAAA\n", encoding="utf-8")
            output.write_text("alignment\n", encoding="utf-8")
            spec = self.make_spec(source)
            snapshot = capture_stage(spec)
            write_stage_manifest(snapshot, manifest, {"paf": output})
            output.write_text("corrupt\n", encoding="utf-8")
            corrupt = validate_stage_cache(snapshot, manifest, {"paf": output})
            output.unlink()
            missing = validate_stage_cache(snapshot, manifest, {"paf": output})

        self.assertFalse(corrupt.hit)
        self.assertIn(corrupt.reason, {"output_size_mismatch", "output_hash_mismatch"})
        self.assertFalse(missing.hit)
        self.assertEqual(missing.reason, "output_missing")

    def test_existing_output_without_valid_manifest_is_not_reused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.fa"
            output = root / "output.paf"
            manifest = root / "manifest.json"
            source.write_text(">u1\nAAAA\n", encoding="utf-8")
            output.write_text("stale\n", encoding="utf-8")
            snapshot = capture_stage(self.make_spec(source))
            missing = validate_stage_cache(snapshot, manifest, {"paf": output})
            manifest.write_text("not-json\n", encoding="utf-8")
            invalid = validate_stage_cache(snapshot, manifest, {"paf": output})

        self.assertEqual(missing.reason, "manifest_missing")
        self.assertEqual(invalid.reason, "manifest_invalid")

    def test_manifest_is_not_written_if_input_changes_during_stage(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.fa"
            output = root / "output.paf"
            manifest = root / "manifest.json"
            source.write_text(">u1\nAAAA\n", encoding="utf-8")
            snapshot = capture_stage(self.make_spec(source))
            output.write_text("alignment\n", encoding="utf-8")
            source.write_text(">u1\nCCCC\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "inputs changed"):
                write_stage_manifest(snapshot, manifest, {"paf": output})
            manifest_exists = manifest.exists()

        self.assertFalse(manifest_exists)

    def test_invalidation_removes_only_completion_marker(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            output = root / "output.paf"
            manifest.write_text("marker\n", encoding="utf-8")
            output.write_text("data\n", encoding="utf-8")
            invalidate_stage_manifest(manifest)
            manifest_exists = manifest.exists()
            output_exists = output.exists()

        self.assertFalse(manifest_exists)
        self.assertTrue(output_exists)

    def test_atomic_writer_preserves_previous_output_on_failure(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "output.paf"
            output.write_text("previous\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "failed"):
                with atomic_text_writer(output) as handle:
                    handle.write("partial\n")
                    raise RuntimeError("failed")
            observed = output.read_text(encoding="utf-8")
            temporary_files = tuple(root.glob(".output.paf.*.tmp"))

        self.assertEqual(observed, "previous\n")
        self.assertEqual(temporary_files, ())


if __name__ == "__main__":
    unittest.main()
