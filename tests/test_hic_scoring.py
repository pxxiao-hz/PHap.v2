from __future__ import annotations

import unittest

from phap_core.hic_scoring import (
    HiCGroupAssignment,
    ScoreMode,
    assign_candidate_groups,
    calculate_group_scores,
    validate_pair_links,
)


class HiCScoringTests(unittest.TestCase):
    def setUp(self) -> None:
        self.groups = {
            "h1": {"anchor_1"},
            "h2": {"anchor_2"},
            "h3": {"anchor_3"},
        }
        self.links = {"anchor_1": 100.0, "anchor_2": 60.0, "anchor_3": 20.0}
        self.re_sites = {"anchor_1": 100, "anchor_2": 10, "anchor_3": 10}

    def assign(
        self,
        *,
        mode: ScoreMode,
        dosage: int = 1,
        min_score: float = 0.0,
        min_margin: float = 0.0,
    ) -> HiCGroupAssignment:
        return assign_candidate_groups(
            "candidate",
            self.groups,
            self.links,
            self.re_sites,
            dosage=dosage,
            ploidy=3,
            mode=mode,
            min_score=min_score,
            min_margin=min_margin,
        )

    def test_raw_and_re_density_can_choose_different_groups(self) -> None:
        raw = self.assign(mode="raw")
        normalized = self.assign(mode="re_density")

        self.assertEqual(raw.selected_groups, ("h1",))
        self.assertEqual(normalized.selected_groups, ("h2",))
        self.assertEqual(raw.top_margin, 40.0)
        self.assertEqual(normalized.top_margin, 4.0)

        audit = {score.group_id: score for score in raw.scores}
        self.assertEqual(audit["h1"].raw_score, 100.0)
        self.assertEqual(audit["h1"].re_density, 1.0)
        self.assertEqual(audit["h1"].raw_rank, 1)
        self.assertEqual(audit["h1"].re_density_rank, 3)
        self.assertEqual(audit["h2"].raw_rank, 2)
        self.assertEqual(audit["h2"].re_density_rank, 1)

    def test_zero_re_denominator_is_unsupported_not_replaced_with_one(self) -> None:
        re_sites = dict(self.re_sites)
        re_sites["anchor_2"] = 0
        scores = calculate_group_scores("candidate", self.groups, self.links, re_sites)
        h2 = next(score for score in scores if score.group_id == "h2")

        self.assertEqual(h2.re_sites, 0)
        self.assertIsNone(h2.re_density)
        self.assertIsNone(h2.re_density_rank)
        self.assertEqual(h2.re_status, "zero_re_sites")

        assignment = assign_candidate_groups(
            "candidate",
            self.groups,
            self.links,
            re_sites,
            dosage=1,
            ploidy=3,
            mode="re_density",
            min_score=0.0,
            min_margin=0.0,
        )
        self.assertEqual(assignment.status, "ambiguous")
        self.assertEqual(assignment.reason, "unsupported_re_density")
        self.assertEqual(assignment.selected_groups, ())

        raw_assignment = assign_candidate_groups(
            "candidate",
            self.groups,
            self.links,
            re_sites,
            dosage=1,
            ploidy=3,
            mode="raw",
            min_score=0.0,
            min_margin=0.0,
        )
        self.assertEqual(raw_assignment.status, "assigned")
        self.assertEqual(raw_assignment.selected_groups, ("h1",))

    def test_dosage_selects_exactly_d_groups(self) -> None:
        assignment = self.assign(mode="raw", dosage=2, min_margin=1.0)

        self.assertEqual(assignment.status, "assigned")
        self.assertEqual(assignment.selected_groups, ("h1", "h2"))
        self.assertEqual(len(assignment.selected_groups), assignment.dosage)
        self.assertEqual(assignment.selection_margin, 40.0)

    def test_tie_at_dosage_boundary_is_ambiguous_despite_stable_tiebreak(self) -> None:
        links = {"anchor_1": 50.0, "anchor_2": 50.0, "anchor_3": 10.0}
        assignment = assign_candidate_groups(
            "candidate",
            dict(reversed(tuple(self.groups.items()))),
            links,
            self.re_sites,
            dosage=1,
            ploidy=3,
            mode="raw",
            min_score=1.0,
            min_margin=0.0,
        )

        self.assertEqual(assignment.status, "ambiguous")
        self.assertEqual(assignment.reason, "non_unique_boundary")
        self.assertEqual(assignment.selected_groups, ())
        self.assertEqual(assignment.top_margin, 0.0)
        ranks = {score.group_id: score.raw_rank for score in assignment.scores}
        self.assertEqual(ranks, {"h1": 1, "h2": 2, "h3": 3})
        self.assertEqual(tuple(score.group_id for score in assignment.scores), ("h1", "h2", "h3"))

    def test_ploidy_is_not_hard_coded_to_four(self) -> None:
        groups = {f"h{index}": {f"a{index}"} for index in range(1, 7)}
        links = {f"a{index}": float(70 - index * 10) for index in range(1, 7)}
        re_sites = {f"a{index}": 10 for index in range(1, 7)}

        assignment = assign_candidate_groups(
            "candidate",
            groups,
            links,
            re_sites,
            dosage=3,
            ploidy=6,
            mode="raw",
            min_score=1.0,
            min_margin=1.0,
        )

        self.assertEqual(assignment.status, "assigned")
        self.assertEqual(assignment.selected_groups, ("h1", "h2", "h3"))
        self.assertEqual(len(assignment.scores), 6)

    def test_low_support_is_ambiguous(self) -> None:
        assignment = self.assign(mode="raw", dosage=2, min_score=70.0)

        self.assertEqual(assignment.status, "ambiguous")
        self.assertEqual(assignment.reason, "below_min_score")
        self.assertEqual(assignment.selected_groups, ())

    def test_margin_threshold_is_explicit_and_applied_at_selection_boundary(self) -> None:
        assignment = self.assign(mode="raw", dosage=1, min_score=1.0, min_margin=50.0)

        self.assertEqual(assignment.status, "ambiguous")
        self.assertEqual(assignment.reason, "below_min_margin")
        self.assertEqual(assignment.selection_margin, 40.0)

    def test_dosage_equal_to_ploidy_requires_no_hic_ranking_guess(self) -> None:
        groups = {"h1": {"a1"}, "h2": {"a2"}}
        assignment = assign_candidate_groups(
            "candidate",
            groups,
            {},
            {"a1": 0, "a2": 0},
            dosage=2,
            ploidy=2,
            mode="re_density",
            min_score=100.0,
            min_margin=100.0,
        )

        self.assertEqual(assignment.status, "assigned")
        self.assertEqual(assignment.reason, "dosage_equals_ploidy")
        self.assertEqual(assignment.selected_groups, ("h1", "h2"))
        self.assertTrue(all(score.re_status == "zero_re_sites" for score in assignment.scores))

    def test_reverse_orientations_cannot_double_count_one_hic_pair(self) -> None:
        with self.assertRaisesRegex(ValueError, "both orientations"):
            validate_pair_links(
                {
                    ("candidate", "anchor_1"): 10,
                    ("anchor_1", "candidate"): 10,
                }
            )


if __name__ == "__main__":
    unittest.main()
