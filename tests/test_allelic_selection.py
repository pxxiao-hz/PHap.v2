from __future__ import annotations

import unittest

from phap_core.allelic_bins import AllelicBinCandidate, AllelicBinKey
from phap_core.allelic_selection import select_dosage_capacity_bins


class AllelicSelectionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.key = AllelicBinKey("chr1", 0, 100)

    def candidate(self, unitig_id: str, support: int) -> AllelicBinCandidate:
        return AllelicBinCandidate(self.key, unitig_id, support)

    def test_capacity_selector_recovers_a_valid_combination(self) -> None:
        decision = select_dosage_capacity_bins(
            [
                self.candidate("d3", 100),
                self.candidate("d2", 99),
                self.candidate("d1", 98),
                self.candidate("x", 97),
            ],
            {"d3": 3, "d2": 2, "d1": 1, "x": 1},
            {"d3": "triplotig", "d2": "diplotig", "d1": "haplotig", "x": "haplotig"},
            ploidy=4,
        )[0]
        self.assertEqual(decision.status, "selected")
        self.assertEqual(decision.selected_unitigs, ("d1", "d3"))
        self.assertEqual(decision.selected_copy_count, 4)
        self.assertEqual(decision.best_objective, 398)

    def test_selector_uses_capacity_dp_instead_of_greedy_ranking(self) -> None:
        decision = select_dosage_capacity_bins(
            [
                self.candidate("d3", 90),
                self.candidate("d2a", 80),
                self.candidate("d2b", 79),
            ],
            {"d3": 3, "d2a": 2, "d2b": 2},
            {},
            ploidy=4,
        )[0]
        self.assertEqual(decision.selected_unitigs, ("d2a", "d2b"))
        self.assertEqual(decision.best_objective, 318)

    def test_equal_best_subsets_are_ambiguous_not_id_tiebroken(self) -> None:
        candidates = [self.candidate("collapsed", 100)] + [
            self.candidate(f"h{index}", 100) for index in range(1, 5)
        ]
        decision = select_dosage_capacity_bins(
            candidates,
            {"collapsed": 4, **{f"h{index}": 1 for index in range(1, 5)}},
            {},
            ploidy=4,
        )[0]
        self.assertEqual(decision.status, "ambiguous")
        self.assertEqual(decision.selected_unitigs, ())
        self.assertEqual(decision.best_objective, 400)
        self.assertEqual(decision.optimal_solution_count, 2)
        self.assertTrue(all(row.status == "ambiguous" for row in decision.candidates))

    def test_invalid_dosages_are_unassigned_and_do_not_consume_capacity(self) -> None:
        decision = select_dosage_capacity_bins(
            [
                self.candidate("valid", 80),
                self.candidate("low", 100),
                self.candidate("high", 100),
                self.candidate("missing", 100),
            ],
            {"valid": 1, "low": None, "high": 5},
            {"valid": "haplotig", "low": "low_coverage", "high": "high_copy"},
            ploidy=4,
        )[0]
        self.assertEqual(decision.status, "partial_capacity")
        self.assertEqual(decision.selected_unitigs, ("valid",))
        reasons = {row.candidate.unitig_id: row.reason for row in decision.candidates}
        self.assertEqual(reasons["low"], "invalid_or_missing_dosage")
        self.assertEqual(reasons["high"], "dosage_outside_ploidy")
        self.assertEqual(reasons["missing"], "invalid_or_missing_dosage")

    def test_optimal_solution_count_is_exact(self) -> None:
        decision = select_dosage_capacity_bins(
            [self.candidate(f"d2_{index}", 50) for index in range(3)],
            {f"d2_{index}": 2 for index in range(3)},
            {},
            ploidy=4,
        )[0]
        self.assertEqual(decision.status, "ambiguous")
        self.assertEqual(decision.optimal_solution_count, 3)

    def test_explicit_bin_support_threshold_routes_slivers_to_unassigned(self) -> None:
        decision = select_dosage_capacity_bins(
            [self.candidate("sliver", 1), self.candidate("supported", 60)],
            {"sliver": 1, "supported": 1},
            {},
            ploidy=2,
            min_support_bases=2,
            min_bin_coverage=0.5,
        )[0]
        self.assertEqual(decision.selected_unitigs, ("supported",))
        reasons = {row.candidate.unitig_id: row.reason for row in decision.candidates}
        self.assertEqual(reasons["sliver"], "insufficient_support_bases")

    def test_ploidy_six_and_candidate_order_are_supported(self) -> None:
        candidates = [
            self.candidate("d3", 100),
            self.candidate("d2", 100),
            self.candidate("d1", 100),
        ]
        dosage = {"d3": 3, "d2": 2, "d1": 1}
        forward = select_dosage_capacity_bins(candidates, dosage, {}, ploidy=6)
        reverse = select_dosage_capacity_bins(
            list(reversed(candidates)), dosage, {}, ploidy=6
        )
        self.assertEqual(forward, reverse)
        self.assertEqual(forward[0].selected_unitigs, ("d1", "d2", "d3"))
        self.assertEqual(forward[0].selected_copy_count, 6)


if __name__ == "__main__":
    unittest.main()
