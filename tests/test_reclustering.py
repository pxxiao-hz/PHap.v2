from __future__ import annotations

import unittest

from phap_core.reclustering import reassign_unitigs


class ReclusterTests(unittest.TestCase):
    def test_assigns_exact_dosage_within_one_family(self) -> None:
        result = reassign_unitigs(
            ["candidate"],
            {
                "a1": {"chr1_group1"},
                "a2": {"chr1_group2"},
                "a3": {"chr1_group3"},
            },
            {"a1": 1, "a2": 1, "a3": 1, "candidate": 2},
            {
                "a1": "haplotig",
                "a2": "haplotig",
                "a3": "haplotig",
                "candidate": "diplotig",
            },
            {
                ("candidate", "a1"): 30,
                ("candidate", "a2"): 20,
                ("candidate", "a3"): 1,
            },
            {"a1": 10, "a2": 10, "a3": 10},
            ploidy=3,
            score_mode="raw",
            min_score=1,
            min_margin=1,
            known_locus="chr1",
        )
        decision = result.decisions[0]
        self.assertEqual(decision.status, "assigned")
        self.assertEqual(
            decision.selected_groups,
            ("chr1_group1", "chr1_group2"),
        )

    def test_no_hic_with_known_locus_is_retained_unresolved(self) -> None:
        result = reassign_unitigs(
            ["candidate"],
            {"a1": {"group1"}, "a2": {"group2"}},
            {"a1": 1, "a2": 1, "candidate": 1},
            {"a1": "haplotig", "a2": "haplotig", "candidate": "haplotig"},
            {},
            {"a1": 10, "a2": 10},
            ploidy=2,
            score_mode="raw",
            min_score=0,
            min_margin=0,
            known_locus="chr1",
        )
        decision = result.decisions[0]
        self.assertEqual(decision.status, "locus_assigned_haplotype_unresolved")
        self.assertEqual(decision.reason, "no_hic_support")
        self.assertEqual(decision.selected_groups, ())

    def test_no_hic_without_locus_is_not_assigned(self) -> None:
        result = reassign_unitigs(
            ["candidate"],
            {"a1": {"chr1_group1"}, "a2": {"chr1_group2"}},
            {"a1": 1, "a2": 1, "candidate": 1},
            {"a1": "haplotig", "a2": "haplotig", "candidate": "haplotig"},
            {},
            {"a1": 10, "a2": 10},
            ploidy=2,
            score_mode="raw",
            min_score=0,
            min_margin=0,
        )
        decision = result.decisions[0]
        self.assertEqual(decision.status, "unassigned")
        self.assertEqual(decision.selected_groups, ())

    def test_competing_chromosome_families_are_ambiguous(self) -> None:
        result = reassign_unitigs(
            ["candidate"],
            {
                "a11": {"chr1_group1"},
                "a12": {"chr1_group2"},
                "a21": {"chr2_group1"},
                "a22": {"chr2_group2"},
            },
            {
                "a11": 1,
                "a12": 1,
                "a21": 1,
                "a22": 1,
                "candidate": 1,
            },
            {
                "a11": "haplotig",
                "a12": "haplotig",
                "a21": "haplotig",
                "a22": "haplotig",
                "candidate": "haplotig",
            },
            {
                ("candidate", "a11"): 10,
                ("candidate", "a21"): 10,
            },
            {"a11": 10, "a12": 10, "a21": 10, "a22": 10},
            ploidy=2,
            score_mode="raw",
            min_score=1,
            min_margin=0,
        )
        decision = result.decisions[0]
        self.assertEqual(decision.status, "ambiguous")
        self.assertEqual(decision.reason, "competing_chromosome_families")

    def test_low_coverage_never_defaults_to_dosage_one(self) -> None:
        result = reassign_unitigs(
            ["low"],
            {"a1": {"group1"}, "a2": {"group2"}},
            {"a1": 1, "a2": 1, "low": None},
            {"a1": "haplotig", "a2": "haplotig", "low": "low_coverage"},
            {},
            {"a1": 10, "a2": 10},
            ploidy=2,
            score_mode="raw",
            min_score=0,
            min_margin=0,
            known_locus="chr1",
        )
        decision = result.decisions[0]
        self.assertEqual(decision.status, "locus_assigned_haplotype_unresolved")
        self.assertEqual(decision.reason, "invalid_or_missing_dosage")

    def test_six_groups_are_supported(self) -> None:
        memberships = {
            f"a{index}": {f"group{index}"} for index in range(1, 7)
        }
        links = {
            ("candidate", f"a{index}"): float(10 - index)
            for index in range(1, 7)
        }
        result = reassign_unitigs(
            ["candidate"],
            memberships,
            {
                **{f"a{index}": 1 for index in range(1, 7)},
                "candidate": 3,
            },
            {
                **{f"a{index}": "haplotig" for index in range(1, 7)},
                "candidate": "dosage_3",
            },
            links,
            {f"a{index}": 10 for index in range(1, 7)},
            ploidy=6,
            score_mode="raw",
            min_score=1,
            min_margin=1,
        )
        self.assertEqual(
            result.decisions[0].selected_groups,
            ("group1", "group2", "group3"),
        )

    def test_preassigned_copy_count_mismatch_stops_the_stage(self) -> None:
        with self.assertRaisesRegex(ValueError, "expected dosage 2"):
            reassign_unitigs(
                ["candidate"],
                {"collapsed": {"group1"}},
                {"collapsed": 2, "candidate": 1},
                {"collapsed": "diplotig", "candidate": "haplotig"},
                {},
                {"collapsed": 10, "candidate": 10},
                ploidy=2,
                score_mode="raw",
                min_score=0,
                min_margin=0,
                declared_groups=("group1", "group2"),
            )


if __name__ == "__main__":
    unittest.main()
