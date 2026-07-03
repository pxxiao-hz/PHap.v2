from __future__ import annotations

import unittest

from phap_core.alignment_chains import (
    select_longest_targets,
    summarize_paf_chains,
)
from phap_core.paf import parse_paf_lines


class AlignmentChainTests(unittest.TestCase):
    def test_overlapping_query_intervals_use_union_coverage(self) -> None:
        records = parse_paf_lines(
            [
                "u\t1000\t0\t600\t+\tchr1\t100000\t10000\t10600"
                "\t570\t600\t60\ttp:A:P\n",
                "u\t1000\t400\t1000\t+\tchr1\t100000\t10400\t11000"
                "\t570\t600\t60\ttp:A:P\n",
            ]
        )
        chain = summarize_paf_chains(
            records,
            max_query_gap=1000,
            max_target_gap=1000,
        )[0]
        self.assertEqual(chain.union_query_bases, 1000)
        self.assertEqual(chain.query_coverage, 1.0)
        self.assertEqual(chain.identity, 0.95)
        self.assertEqual(chain.score, 0.95)
        self.assertTrue(chain.collinear)

    def test_reverse_chain_requires_decreasing_target_coordinates(self) -> None:
        valid = parse_paf_lines(
            [
                "u\t1000\t0\t400\t-\tchr3\t20000\t9600\t10000"
                "\t380\t400\t60\ttp:A:P\n",
                "u\t1000\t400\t900\t-\tchr3\t20000\t9100\t9600"
                "\t475\t500\t60\ttp:A:P\n",
            ]
        )
        invalid = parse_paf_lines(
            [
                "bad\t1000\t0\t400\t-\tchr3\t20000\t9600\t10000"
                "\t380\t400\t60\ttp:A:P\n",
                "bad\t1000\t400\t900\t-\tchr3\t20000\t10000\t10500"
                "\t475\t500\t60\ttp:A:P\n",
            ]
        )
        valid_chain = summarize_paf_chains(
            valid,
            max_query_gap=1000,
            max_target_gap=1000,
        )[0]
        invalid_chain = summarize_paf_chains(
            invalid,
            max_query_gap=1000,
            max_target_gap=1000,
        )[0]
        self.assertTrue(valid_chain.collinear)
        self.assertEqual(valid_chain.query_coverage, 0.9)
        self.assertEqual(valid_chain.identity, 0.95)
        self.assertFalse(invalid_chain.collinear)

    def test_mixed_orientation_is_explicitly_non_collinear(self) -> None:
        records = parse_paf_lines(
            [
                "u\t1000\t0\t400\t+\tchr1\t10000\t0\t400"
                "\t400\t400\t60\ttp:A:P\n",
                "u\t1000\t400\t800\t-\tchr1\t10000\t400\t800"
                "\t400\t400\t60\ttp:A:P\n",
            ]
        )
        chain = summarize_paf_chains(
            records,
            max_query_gap=1000,
            max_target_gap=1000,
        )[0]
        self.assertEqual(chain.strand, "mixed")
        self.assertFalse(chain.collinear)

    def test_gap_thresholds_are_applied_in_bases(self) -> None:
        records = parse_paf_lines(
            [
                "u\t2000\t0\t400\t+\tchr1\t10000\t0\t400"
                "\t400\t400\t60\ttp:A:P\n",
                "u\t2000\t1400\t1800\t+\tchr1\t10000\t1400\t1800"
                "\t400\t400\t60\ttp:A:P\n",
            ]
        )
        chain = summarize_paf_chains(
            records,
            max_query_gap=999,
            max_target_gap=999,
        )[0]
        self.assertFalse(chain.collinear)

    def test_target_selection_is_length_then_id_deterministic(self) -> None:
        records = parse_paf_lines(
            [
                "u1\t100\t0\t80\t+\tchrB\t1000\t0\t80\t80\t80\t60\ttp:A:P\n",
                "u2\t100\t0\t80\t+\tchrA\t1000\t0\t80\t80\t80\t60\ttp:A:P\n",
                "u3\t100\t0\t80\t+\tshort\t500\t0\t80\t80\t80\t60\ttp:A:P\n",
            ]
        )
        self.assertEqual(
            select_longest_targets(records, target_count=2),
            ("chrA", "chrB"),
        )


if __name__ == "__main__":
    unittest.main()
