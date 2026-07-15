from __future__ import annotations

import unittest

from phap_core.allelic_bins import (
    AllelicBinCandidate,
    AllelicBinKey,
    TargetIntervalEvidence,
    build_allelic_bin_candidates,
    select_legacy_top_n_candidates,
)


class AllelicBinConstructionTests(unittest.TestCase):
    def test_overlapping_alignments_use_interval_union_per_bin(self) -> None:
        candidates = build_allelic_bin_candidates(
            [
                TargetIntervalEvidence("u1", "chr1", 100, 0, 80),
                TargetIntervalEvidence("u1", "chr1", 100, 20, 100),
            ],
            bin_size=100,
        )
        self.assertEqual(
            candidates,
            (AllelicBinCandidate(AllelicBinKey("chr1", 0, 100), "u1", 100),),
        )

    def test_alignment_is_clipped_into_half_open_target_bins(self) -> None:
        candidates = build_allelic_bin_candidates(
            [TargetIntervalEvidence("u1", "chr1", 250, 80, 230)],
            bin_size=100,
        )
        self.assertEqual(
            candidates,
            (
                AllelicBinCandidate(AllelicBinKey("chr1", 0, 100), "u1", 20),
                AllelicBinCandidate(AllelicBinKey("chr1", 100, 200), "u1", 100),
                AllelicBinCandidate(AllelicBinKey("chr1", 200, 250), "u1", 30),
            ),
        )

    def test_legacy_top_n_is_support_ranked_with_stable_id_ties(self) -> None:
        key = AllelicBinKey("chr1", 0, 100)
        selected = select_legacy_top_n_candidates(
            [
                AllelicBinCandidate(key, "u3", 80),
                AllelicBinCandidate(key, "u2", 90),
                AllelicBinCandidate(key, "u1", 90),
            ],
            top_n=2,
        )
        self.assertEqual(
            selected,
            (
                AllelicBinCandidate(key, "u1", 90),
                AllelicBinCandidate(key, "u2", 90),
            ),
        )

    def test_candidate_output_is_independent_of_alignment_order(self) -> None:
        alignments = [
            TargetIntervalEvidence("u2", "chr2", 100, 0, 50),
            TargetIntervalEvidence("u1", "chr1", 100, 20, 100),
            TargetIntervalEvidence("u1", "chr1", 100, 0, 40),
        ]
        forward = build_allelic_bin_candidates(alignments, bin_size=100)
        reverse = build_allelic_bin_candidates(reversed(alignments), bin_size=100)
        self.assertEqual(forward, reverse)


if __name__ == "__main__":
    unittest.main()
