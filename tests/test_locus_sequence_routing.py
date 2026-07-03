from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from phap_core.fasta import FastaRecord, read_fasta
from utils.extract_locus_sequences import partition_records, write_partition


class LocusSequenceRoutingTests(unittest.TestCase):
    def test_partition_is_complete_disjoint_and_preserves_headers(self) -> None:
        records = (
            FastaRecord("u1", "u1 description", "AAAA"),
            FastaRecord("u2", "u2", "CCCC"),
            FastaRecord("u3", "u3", "GGGG"),
        )
        decisions = {
            "u1": (
                "locus_assigned_haplotype_unresolved",
                "chr1",
                "mT2T_supports_locus_not_haplotype",
            ),
            "u2": ("ambiguous", "", "competing_mT2T_loci"),
            "u3": ("unassigned", "", "no_collinear_locus_meets_thresholds"),
        }
        by_locus, unplaced, routing = partition_records(records, decisions)
        self.assertEqual(tuple(by_locus), ("chr1",))
        self.assertEqual(by_locus["chr1"][0].header, "u1 description")
        self.assertEqual(
            tuple(record.identifier for record in unplaced),
            ("u2", "u3"),
        )
        self.assertEqual(len(routing), 3)

    def test_outputs_use_original_ids_and_one_destination_each(self) -> None:
        records = (
            FastaRecord("u1", "u1 description", "AAAA"),
            FastaRecord("u2", "u2", "CCCC"),
        )
        decisions = {
            "u1": (
                "locus_assigned_haplotype_unresolved",
                "scaffold_A",
                "mT2T_supports_locus_not_haplotype",
            ),
            "u2": ("ambiguous", "", "competing_mT2T_loci"),
        }
        with tempfile.TemporaryDirectory() as directory:
            write_partition(records, decisions, directory)
            assigned = (Path(directory) / "scaffold_A.putg.fa").read_text(
                encoding="utf-8"
            )
            unplaced = (Path(directory) / "un_chr.fa").read_text(encoding="utf-8")
            routing = (Path(directory) / "locus_sequence_routing.tsv").read_text(
                encoding="utf-8"
            )
        self.assertEqual(assigned, ">u1 description\nAAAA\n")
        self.assertEqual(unplaced, ">u2\nCCCC\n")
        self.assertIn("u1\tscaffold_A.putg.fa", routing)
        self.assertIn("u2\tun_chr.fa", routing)

    def test_fasta_parser_rejects_duplicate_first_token_ids(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "duplicate.fa"
            fasta.write_text(
                ">u description one\nAAAA\n>u description two\nCCCC\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "duplicate FASTA identifier"):
                read_fasta(fasta)

    def test_rerun_removes_only_manifest_owned_stale_locus_files(self) -> None:
        records = (FastaRecord("u1", "u1", "AAAA"),)
        with tempfile.TemporaryDirectory() as directory:
            write_partition(
                records,
                {"u1": ("assigned", "chrA", "first")},
                directory,
            )
            write_partition(
                records,
                {"u1": ("assigned", "chrB", "second")},
                directory,
            )
            old_exists = (Path(directory) / "chrA.putg.fa").exists()
            new_exists = (Path(directory) / "chrB.putg.fa").exists()
        self.assertFalse(old_exists)
        self.assertTrue(new_exists)


if __name__ == "__main__":
    unittest.main()
