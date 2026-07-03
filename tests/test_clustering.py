from __future__ import annotations

import unittest

from phap_core.clustering import AllelicBin, cluster_allelic_bins


class ClusteringTests(unittest.TestCase):
    def test_first_valid_bin_seeds_all_anonymous_groups_deterministically(self) -> None:
        result = cluster_allelic_bins(
            [AllelicBin("chr1", 0, 100, ("d", "b", "a", "c"))],
            {"a": 1, "b": 1, "c": 1, "d": 1},
            {"a": "haplotig", "b": "haplotig", "c": "haplotig", "d": "haplotig"},
            {
                ("a", "external_a"): 1,
                ("b", "external_b"): 1,
                ("c", "external_c"): 1,
                ("d", "external_d"): 1,
            },
            {"a": 1, "b": 1, "c": 1, "d": 1},
            ploidy=4,
            score_mode="raw",
            min_score=0,
            min_margin=0,
        )
        self.assertEqual(
            dict(result.group_unitigs),
            {
                "group1": ("a",),
                "group2": ("b",),
                "group3": ("c",),
                "group4": ("d",),
            },
        )
        self.assertTrue(
            all(
                decision.reason == "deterministic_initial_seed"
                for decision in result.decisions
            )
        )

    def test_dynamic_ploidy_and_dosage_are_conserved(self) -> None:
        rows = [
            AllelicBin("chr1", 0, 100, ("collapsed", "hap")),
            AllelicBin("chr1", 100, 200, ("collapsed", "hap")),
        ]
        result = cluster_allelic_bins(
            rows,
            {"collapsed": 2, "hap": 1},
            {"collapsed": "diplotig", "hap": "haplotig"},
            {
                ("collapsed", "external_collapsed"): 1,
                ("hap", "external_hap"): 1,
            },
            {"collapsed": 2, "hap": 2},
            ploidy=3,
            score_mode="raw",
            min_score=0,
            min_margin=0,
        )
        groups = dict(result.group_unitigs)
        self.assertEqual(groups["group1"], ("collapsed",))
        self.assertEqual(groups["group2"], ("collapsed",))
        self.assertEqual(groups["group3"], ("hap",))

    def test_raw_and_density_modes_drive_different_assignments(self) -> None:
        rows = [
            AllelicBin("chr1", 0, 100, ("a1", "a2")),
            AllelicBin("chr1", 100, 200, ("candidate",)),
        ]
        dosage = {"a1": 1, "a2": 1, "candidate": 1}
        states = {"a1": "haplotig", "a2": "haplotig", "candidate": "haplotig"}
        links = {("candidate", "a1"): 100.0, ("candidate", "a2"): 60.0}
        links.update({("a1", "external_1"): 1.0, ("a2", "external_2"): 1.0})
        re_sites = {"a1": 100, "a2": 10, "candidate": 10}
        raw = cluster_allelic_bins(
            rows,
            dosage,
            states,
            links,
            re_sites,
            ploidy=2,
            score_mode="raw",
            min_score=0,
            min_margin=0,
        )
        density = cluster_allelic_bins(
            rows,
            dosage,
            states,
            links,
            re_sites,
            ploidy=2,
            score_mode="re_density",
            min_score=0,
            min_margin=0,
        )
        raw_decision = next(
            decision for decision in raw.decisions if decision.unitig_id == "candidate"
        )
        density_decision = next(
            decision
            for decision in density.decisions
            if decision.unitig_id == "candidate"
        )
        self.assertEqual(raw_decision.selected_groups, ("group1",))
        self.assertEqual(density_decision.selected_groups, ("group2",))
        self.assertTrue(raw_decision.scores)
        self.assertTrue(density_decision.scores)

    def test_no_hic_is_locus_assigned_but_not_given_a_group(self) -> None:
        rows = [
            AllelicBin("chr1", 0, 100, ("anchor1", "anchor2")),
            AllelicBin("chr1", 100, 200, ("unsupported",)),
        ]
        result = cluster_allelic_bins(
            rows,
            {"anchor1": 1, "anchor2": 1, "unsupported": 1},
            {"anchor1": "haplotig", "anchor2": "haplotig", "unsupported": "haplotig"},
            {
                ("anchor1", "external_1"): 1,
                ("anchor2", "external_2"): 1,
            },
            {"anchor1": 10, "anchor2": 10, "unsupported": 10},
            ploidy=2,
            score_mode="raw",
            min_score=0,
            min_margin=0,
        )
        decision = next(
            decision
            for decision in result.decisions
            if decision.unitig_id == "unsupported"
        )
        self.assertEqual(decision.status, "locus_assigned_haplotype_unresolved")
        self.assertEqual(decision.reason, "no_hic_support")
        self.assertEqual(decision.selected_groups, ())
        self.assertNotIn(
            "unsupported",
            {unitig for _, unitigs in result.group_unitigs for unitig in unitigs},
        )

    def test_low_coverage_is_not_silently_treated_as_haplotig(self) -> None:
        rows = [AllelicBin("chr1", 0, 100, ("low",))]
        result = cluster_allelic_bins(
            rows,
            {"low": None},
            {"low": "low_coverage"},
            {},
            {"low": 1},
            ploidy=4,
            score_mode="raw",
            min_score=0,
            min_margin=0,
        )
        decision = result.decisions[0]
        self.assertEqual(decision.source_state, "low_coverage")
        self.assertEqual(decision.status, "locus_assigned_haplotype_unresolved")
        self.assertEqual(decision.reason, "invalid_or_missing_dosage")
        self.assertEqual(decision.selected_groups, ())

    def test_insufficient_capacity_never_creates_partial_copy_assignment(self) -> None:
        rows = [
            AllelicBin("chr1", 0, 100, ("a1", "a2", "a3", "a4")),
            AllelicBin("chr1", 100, 200, ("u2", "u3")),
        ]
        result = cluster_allelic_bins(
            rows,
            {"a1": 1, "a2": 1, "a3": 1, "a4": 1, "u2": 2, "u3": 3},
            {
                "a1": "haplotig",
                "a2": "haplotig",
                "a3": "haplotig",
                "a4": "haplotig",
                "u2": "diplotig",
                "u3": "triplotig",
            },
            {
                ("a1", "external_1"): 1,
                ("a2", "external_2"): 1,
                ("a3", "external_3"): 1,
                ("a4", "external_4"): 1,
                ("u3", "a1"): 10,
                ("u3", "a2"): 9,
                ("u3", "a3"): 8,
                ("u3", "a4"): 1,
            },
            {"a1": 2, "a2": 2, "a3": 2, "a4": 2, "u2": 2, "u3": 2},
            ploidy=4,
            score_mode="raw",
            min_score=0,
            min_margin=0,
        )
        decisions = {decision.unitig_id: decision for decision in result.decisions}
        self.assertEqual(len(decisions["u3"].selected_groups), 3)
        self.assertEqual(decisions["u2"].selected_groups, ())
        self.assertEqual(decisions["u2"].status, "ambiguous")
        self.assertEqual(decisions["u2"].reason, "insufficient_bin_capacity")

    def test_tied_hic_evidence_is_ambiguous(self) -> None:
        rows = [
            AllelicBin("chr1", 0, 100, ("a1", "a2")),
            AllelicBin("chr1", 100, 200, ("candidate",)),
        ]
        result = cluster_allelic_bins(
            rows,
            {"a1": 1, "a2": 1, "candidate": 1},
            {"a1": "haplotig", "a2": "haplotig", "candidate": "haplotig"},
            {
                ("a1", "external_1"): 1,
                ("a2", "external_2"): 1,
                ("candidate", "a1"): 10,
                ("candidate", "a2"): 10,
            },
            {"a1": 10, "a2": 10, "candidate": 10},
            ploidy=2,
            score_mode="raw",
            min_score=0,
            min_margin=0,
        )
        decision = next(
            decision for decision in result.decisions if decision.unitig_id == "candidate"
        )
        self.assertEqual(decision.status, "ambiguous")
        self.assertEqual(decision.reason, "non_unique_boundary")
        self.assertEqual(decision.selected_groups, ())

    def test_no_hic_unitigs_cannot_seed_haplotype_groups(self) -> None:
        result = cluster_allelic_bins(
            [AllelicBin("chr1", 0, 100, ("a", "b"))],
            {"a": 1, "b": 1},
            {"a": "haplotig", "b": "haplotig"},
            {},
            {"a": 1, "b": 1},
            ploidy=2,
            score_mode="raw",
            min_score=0,
            min_margin=0,
        )
        self.assertEqual(
            {decision.reason for decision in result.decisions},
            {"no_hic_seed_support"},
        )
        self.assertTrue(
            all(
                decision.status == "locus_assigned_haplotype_unresolved"
                for decision in result.decisions
            )
        )
        self.assertTrue(
            all(not unitigs for _, unitigs in result.group_unitigs)
        )


if __name__ == "__main__":
    unittest.main()
