from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from phap_core.recluster_merge import merge_current_recluster_outputs


class ReclusterMergeTests(unittest.TestCase):
    def write_source(self, root: Path, locus_id: str, text: str) -> Path:
        source = root / locus_id / "group.reassignment.cluster.txt"
        source.parent.mkdir(parents=True)
        source.write_text(text, encoding="utf-8")
        return source

    def test_only_current_loci_are_merged_in_stable_order(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "04.recluster"
            self.write_source(root, "old_chr", "old_group\told_u\n")
            self.write_source(root, "chr2", "chr2_group\tu2\n")
            self.write_source(root, "chr1", "chr1_group\tu1\n")
            destination = Path(directory) / "merged.txt"
            audit = Path(directory) / "merge_sources.tsv"

            sources = merge_current_recluster_outputs(
                root,
                ("chr2", "chr1"),
                destination,
                audit,
            )
            merged = destination.read_text(encoding="utf-8")
            audit_text = audit.read_text(encoding="utf-8")

        self.assertEqual(tuple(row.locus_id for row in sources), ("chr1", "chr2"))
        self.assertEqual(merged, "chr1_group\tu1\nchr2_group\tu2\n")
        self.assertNotIn("old_group", merged)
        self.assertIn("chr1/group.reassignment.cluster.txt", audit_text)
        self.assertIn("chr2/group.reassignment.cluster.txt", audit_text)
        self.assertNotIn("old_chr", audit_text)

    def test_empty_current_locus_set_overwrites_stale_merge_with_empty(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "04.recluster"
            self.write_source(root, "old_chr", "old_group\told_u\n")
            destination = Path(directory) / "merged.txt"
            audit = Path(directory) / "merge_sources.tsv"
            destination.write_text("stale\n", encoding="utf-8")

            merge_current_recluster_outputs(root, (), destination, audit)

            self.assertEqual(destination.read_text(encoding="utf-8"), "")
            self.assertEqual(
                audit.read_text(encoding="utf-8"),
                "locus_ID\tsource_file\tsize\tsha256\tdestination\treason\n",
            )

    def test_missing_current_output_preserves_previous_merge(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "04.recluster"
            root.mkdir(parents=True)
            destination = Path(directory) / "merged.txt"
            audit = Path(directory) / "merge_sources.tsv"
            destination.write_text("previous\n", encoding="utf-8")
            audit.write_text("previous audit\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "output is missing"):
                merge_current_recluster_outputs(
                    root,
                    ("chr1",),
                    destination,
                    audit,
                )

            self.assertEqual(destination.read_text(encoding="utf-8"), "previous\n")
            self.assertEqual(audit.read_text(encoding="utf-8"), "previous audit\n")

    def test_duplicate_or_unsafe_locus_is_rejected_before_writes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "04.recluster"
            root.mkdir(parents=True)
            destination = Path(directory) / "merged.txt"
            audit = Path(directory) / "merge_sources.tsv"
            destination.write_text("previous\n", encoding="utf-8")

            for loci in (("chr1", "chr1"), ("../chr1",)):
                with self.subTest(loci=loci):
                    with self.assertRaises(ValueError):
                        merge_current_recluster_outputs(
                            root,
                            loci,
                            destination,
                            audit,
                        )

            self.assertEqual(destination.read_text(encoding="utf-8"), "previous\n")

    def test_mislabeled_group_is_rejected_before_replacing_merge(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "04.recluster"
            self.write_source(root, "chr1", "chr2_group1\tu1\n")
            destination = Path(directory) / "merged.txt"
            audit = Path(directory) / "merge_sources.tsv"
            destination.write_text("previous\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "does not belong to locus"):
                merge_current_recluster_outputs(
                    root,
                    ("chr1",),
                    destination,
                    audit,
                )

            self.assertEqual(destination.read_text(encoding="utf-8"), "previous\n")


if __name__ == "__main__":
    unittest.main()
