from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from phap_core.fasta import FastaRecord
from phap_core.fasta import count_fasta_records
from phap_core.locus_sequence_manifest import (
    locus_id_from_fasta_name,
    read_current_locus_outputs,
)
from utils.extract_locus_sequences import write_partition


class LocusSequenceManifestTests(unittest.TestCase):
    def test_manifest_is_the_sorted_authority_for_current_loci(self) -> None:
        records = (
            FastaRecord("u1", "u1", "AAAA"),
            FastaRecord("u2", "u2", "CCCC"),
            FastaRecord("u3", "u3", "GGGG"),
        )
        decisions = {
            "u1": ("assigned", "chr2", "current"),
            "u2": ("assigned", "chr1", "current"),
            "u3": ("unassigned", "", "no_locus"),
        }
        with tempfile.TemporaryDirectory() as directory:
            write_partition(records, decisions, directory)
            outputs = read_current_locus_outputs(directory)

        self.assertEqual(tuple(row.locus_id for row in outputs), ("chr1", "chr2"))
        self.assertEqual(tuple(row.record_count for row in outputs), (1, 1))

    def test_unmanifested_locus_fasta_is_rejected(self) -> None:
        records = (FastaRecord("u1", "u1", "AAAA"),)
        with tempfile.TemporaryDirectory() as directory:
            write_partition(
                records,
                {"u1": ("assigned", "chr1", "current")},
                directory,
            )
            (Path(directory) / "ghost.putg.fa").write_text(
                ">ghost\nAAAA\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "does not match manifest"):
                read_current_locus_outputs(directory)

    def test_record_count_mismatch_is_rejected(self) -> None:
        records = (FastaRecord("u1", "u1", "AAAA"),)
        with tempfile.TemporaryDirectory() as directory:
            write_partition(
                records,
                {"u1": ("assigned", "chr1", "current")},
                directory,
            )
            manifest = Path(directory) / "locus_sequence_manifest.tsv"
            manifest.write_text(
                "output_file\trecord_count\nchr1.putg.fa\t2\nun_chr.fa\t0\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "record count mismatch"):
                read_current_locus_outputs(directory)

    def test_only_the_terminal_fasta_suffix_is_removed(self) -> None:
        self.assertEqual(
            locus_id_from_fasta_name("chr.putg.fa.copy.putg.fa"),
            "chr.putg.fa.copy",
        )

    def test_streaming_fasta_count_rejects_duplicate_ids(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "duplicate.fa"
            fasta.write_text(
                ">u description one\nAAAA\n>u description two\nCCCC\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "duplicate FASTA identifier"):
                count_fasta_records(fasta)


if __name__ == "__main__":
    unittest.main()
