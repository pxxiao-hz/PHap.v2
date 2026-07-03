from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

from phap_core.reclustering import reassign_unitigs
from utils.recluster_common import (
    parse_group_file,
    read_fasta_with_re_sites,
    write_recluster_outputs,
)


class ReclusterCommonTests(unittest.TestCase):
    def test_both_group_file_stages_have_explicit_parsers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            counted = directory / "counted.groups.tsv"
            counted.write_text(
                "group1\t2\tu1\tu2\n"
                "group2\t0\n",
                encoding="utf-8",
            )
            memberships, groups = parse_group_file(
                str(counted),
                has_count_column=True,
            )
            self.assertEqual(groups, ("group1", "group2"))
            self.assertEqual(memberships, {"u1": {"group1"}, "u2": {"group1"}})

            uncounted = directory / "uncounted.groups.tsv"
            uncounted.write_text(
                "chr1_group1\tu1\tu2\n"
                "chr1_group2\n",
                encoding="utf-8",
            )
            memberships, groups = parse_group_file(
                str(uncounted),
                has_count_column=False,
            )
            self.assertEqual(groups, ("chr1_group1", "chr1_group2"))
            self.assertEqual(
                memberships,
                {"u1": {"chr1_group1"}, "u2": {"chr1_group1"}},
            )

    def test_duplicate_unitig_within_group_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            groups = Path(temporary_directory) / "groups.tsv"
            groups.write_text("group1\tu1\tu1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate unitig"):
                parse_group_file(str(groups), has_count_column=False)

    def test_re_site_counts_do_not_add_a_hidden_pseudocount(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            fasta = Path(temporary_directory) / "unitigs.fa"
            fasta.write_text(
                ">zero\nAAAA\n>two\nGATCGATC\n",
                encoding="utf-8",
            )
            _, lengths, re_sites = read_fasta_with_re_sites(str(fasta), "GATC")
        self.assertEqual(lengths, {"zero": 4, "two": 8})
        self.assertEqual(re_sites, {"zero": 0, "two": 2})

    def test_outputs_keep_unresolved_unitigs_in_the_audit(self) -> None:
        result = reassign_unitigs(
            ["unsupported"],
            {"anchor1": {"group1"}, "anchor2": {"group2"}},
            {"anchor1": 1, "anchor2": 1, "unsupported": 1},
            {
                "anchor1": "haplotig",
                "anchor2": "haplotig",
                "unsupported": "haplotig",
            },
            {},
            {"anchor1": 1, "anchor2": 1, "unsupported": 1},
            ploidy=2,
            score_mode="raw",
            min_score=0,
            min_margin=0,
            known_locus="chr1",
            declared_groups=("group1", "group2"),
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            previous = Path.cwd()
            os.chdir(temporary_directory)
            try:
                write_recluster_outputs(
                    result,
                    {"anchor1": "AAAA", "anchor2": "CCCC", "unsupported": "GGGG"},
                    {"anchor1": 4, "anchor2": 4, "unsupported": 4},
                    {"anchor1": 1, "anchor2": 1, "unsupported": 1},
                    ("group1", "group2"),
                    output_group_id=lambda group: f"chr1_{group}",
                )
                decisions = Path("recluster_decisions.tsv").read_text(
                    encoding="utf-8"
                )
                groups = Path("group.reassignment.cluster.txt").read_text(
                    encoding="utf-8"
                )
                unresolved_fasta = Path("unassigned_unitigs.fa").read_text(
                    encoding="utf-8"
                )
            finally:
                os.chdir(previous)
        self.assertIn(
            "unsupported\thaplotig\t1\tlocus_assigned_haplotype_unresolved",
            decisions,
        )
        self.assertEqual(
            groups,
            "chr1_group1\tanchor1\nchr1_group2\tanchor2\n",
        )
        self.assertEqual(unresolved_fasta, ">unsupported\nGGGG\n")


if __name__ == "__main__":
    unittest.main()
