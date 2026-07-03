from __future__ import annotations

import unittest

from phap_core.locus_rescue import (
    AlignmentChainSummary,
    GroupSpecificEvidence,
    LocusRescueThresholds,
    decide_locus_rescue,
    interval_union_length,
    summarize_locus_candidates,
)


class LocusRescueTests(unittest.TestCase):
    def setUp(self) -> None:
        self.thresholds = LocusRescueThresholds(
            min_identity=0.90,
            min_query_coverage=0.75,
            min_next_best_margin=0.10,
            min_low_coverage_read_support=2,
            min_group_read_support=2,
            min_group_overlap_support=100.0,
        )

    @staticmethod
    def alignment(
        locus_id: str,
        *,
        identity: float = 0.98,
        intervals: tuple[tuple[int, int], ...] = ((0, 900),),
        collinear: bool = True,
    ) -> AlignmentChainSummary:
        return AlignmentChainSummary(
            unitig_id="utg1",
            locus_id=locus_id,
            query_length=1000,
            query_intervals=intervals,
            identity=identity,
            collinear=collinear,
        )

    def test_unique_mT2T_alignment_assigns_locus_but_not_haplotype(self) -> None:
        decision = decide_locus_rescue(
            "utg1",
            [self.alignment("chr01:0-100000")],
            source_state="no_hic_signal",
            independent_read_support=0,
            thresholds=self.thresholds,
        )
        self.assertEqual(decision.status, "locus_assigned_haplotype_unresolved")
        self.assertEqual(decision.assigned_locus, "chr01:0-100000")
        self.assertIsNone(decision.assigned_group)
        self.assertEqual(decision.reason, "mT2T_supports_locus_not_haplotype")

    def test_competing_loci_are_ambiguous(self) -> None:
        decision = decide_locus_rescue(
            "utg1",
            [
                self.alignment("chr01:0-100000", identity=0.98),
                self.alignment("chr02:0-100000", identity=0.97),
            ],
            source_state="no_hic_signal",
            independent_read_support=0,
            thresholds=self.thresholds,
        )
        self.assertEqual(decision.status, "ambiguous")
        self.assertEqual(decision.reason, "competing_mT2T_loci")
        self.assertIsNone(decision.assigned_locus)
        self.assertEqual(decision.best_locus, "chr01:0-100000")
        self.assertEqual(decision.next_best_locus, "chr02:0-100000")

    def test_overlapping_alignment_intervals_use_union_coverage(self) -> None:
        alignment = self.alignment(
            "chr01:0-100000",
            intervals=((0, 600), (400, 1000)),
        )
        candidates = summarize_locus_candidates("utg1", [alignment], self.thresholds)
        self.assertEqual(interval_union_length(alignment.query_intervals), 1000)
        self.assertEqual(candidates[0].union_query_bases, 1000)
        self.assertEqual(candidates[0].query_coverage, 1.0)

    def test_low_coverage_requires_unique_locus_and_independent_reads(self) -> None:
        without_reads = decide_locus_rescue(
            "utg1",
            [self.alignment("chr01:0-100000")],
            source_state="low_coverage",
            independent_read_support=1,
            thresholds=self.thresholds,
        )
        self.assertEqual(without_reads.status, "unassigned")
        self.assertEqual(
            without_reads.reason,
            "low_coverage_without_independent_read_support",
        )

        rescued = decide_locus_rescue(
            "utg1",
            [self.alignment("chr01:0-100000")],
            source_state="low_coverage",
            independent_read_support=2,
            thresholds=self.thresholds,
        )
        self.assertEqual(rescued.status, "locus_assigned_haplotype_unresolved")
        self.assertEqual(rescued.rescue_class, "rescued_haplotig")
        self.assertIsNone(rescued.assigned_group)

    def test_only_extra_group_specific_evidence_can_rescue_group(self) -> None:
        no_group_evidence = decide_locus_rescue(
            "utg1",
            [self.alignment("chr01:0-100000")],
            source_state="no_hic_signal",
            independent_read_support=0,
            thresholds=self.thresholds,
        )
        self.assertEqual(
            no_group_evidence.status,
            "locus_assigned_haplotype_unresolved",
        )

        by_reads = decide_locus_rescue(
            "utg1",
            [self.alignment("chr01:0-100000")],
            source_state="no_hic_signal",
            independent_read_support=0,
            thresholds=self.thresholds,
            group_evidence=[
                GroupSpecificEvidence(
                    group_id="chr01_h2",
                    locus_id="chr01:0-100000",
                    read_support=2,
                )
            ],
        )
        self.assertEqual(by_reads.status, "rescued_group")
        self.assertEqual(by_reads.assigned_group, "chr01_h2")

        by_overlap = decide_locus_rescue(
            "utg1",
            [self.alignment("chr01:0-100000")],
            source_state="no_hic_signal",
            independent_read_support=0,
            thresholds=self.thresholds,
            group_evidence=[
                GroupSpecificEvidence(
                    group_id="chr01_h3",
                    locus_id="chr01:0-100000",
                    overlap_support=100.0,
                )
            ],
        )
        self.assertEqual(by_overlap.status, "rescued_group")
        self.assertEqual(by_overlap.assigned_group, "chr01_h3")

    def test_identity_coverage_and_margin_thresholds_are_inclusive(self) -> None:
        thresholds = LocusRescueThresholds(
            min_identity=0.80,
            min_query_coverage=0.75,
            min_next_best_margin=0.25,
        )
        decision = decide_locus_rescue(
            "utg1",
            [
                self.alignment(
                    "chr01:0-100000",
                    identity=0.80,
                    intervals=((0, 750),),
                ),
                self.alignment(
                    "chr02:0-100000",
                    identity=0.70,
                    intervals=((0, 500),),
                ),
            ],
            source_state="no_hic_signal",
            independent_read_support=0,
            thresholds=thresholds,
        )
        self.assertEqual(decision.status, "locus_assigned_haplotype_unresolved")
        self.assertAlmostEqual(decision.best_identity or 0.0, 0.80)
        self.assertAlmostEqual(decision.best_query_coverage or 0.0, 0.75)
        self.assertAlmostEqual(decision.next_best_margin or 0.0, 0.25)


if __name__ == "__main__":
    unittest.main()
