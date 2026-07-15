from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from phap_core.cluster_tables import write_chromosome_allelic_table


class ClusterTableTests(unittest.TestCase):
    def test_existing_table_is_always_replaced_from_current_source(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "corrected.tsv"
            destination = root / "chr1.tsv"
            source.write_text(
                "chr1\t0\t100\tnew_u1\nchr10\t0\t100\tother\n",
                encoding="utf-8",
            )
            destination.write_text("chr1\t0\t100\tstale\n", encoding="utf-8")
            count = write_chromosome_allelic_table(source, "chr1", destination)
            observed = destination.read_text(encoding="utf-8")

        self.assertEqual(count, 1)
        self.assertEqual(observed, "chr1\t0\t100\tnew_u1\n")

    def test_absent_chromosome_atomically_replaces_stale_table_with_empty(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "corrected.tsv"
            destination = root / "chr1.tsv"
            source.write_text("chr2\t0\t100\tu2\n", encoding="utf-8")
            destination.write_text("chr1\t0\t100\tstale\n", encoding="utf-8")
            count = write_chromosome_allelic_table(source, "chr1", destination)
            observed = destination.read_text(encoding="utf-8")

        self.assertEqual(count, 0)
        self.assertEqual(observed, "")

    def test_malformed_source_does_not_replace_previous_destination(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "corrected.tsv"
            destination = root / "chr1.tsv"
            source.write_text("chr1\t0\t100\n", encoding="utf-8")
            destination.write_text("previous\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "expected chromosome"):
                write_chromosome_allelic_table(source, "chr1", destination)
            observed = destination.read_text(encoding="utf-8")

        self.assertEqual(observed, "previous\n")


if __name__ == "__main__":
    unittest.main()
