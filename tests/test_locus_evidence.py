from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from phap_core.alignment_chains import select_longest_targets
from phap_core.locus_evidence import evaluate_locus_evidence
from phap_core.locus_rescue import LocusRescueThresholds
from phap_core.paf import parse_paf_lines
from utils.paf_locus_evidence import write_outputs


class LocusEvidenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.thresholds = LocusRescueThresholds(
            min_identity=0.8,
            min_query_coverage=0.5,
            min_next_best_margin=0.03,
            min_low_coverage_read_support=2,
        )

    def test_competing_loci_use_complete_paf_and_union_coverage(self) -> None:
        records = parse_paf_lines(
            [
                "u\t1000\t0\t800\t+\tchr1\t100000\t1000\t1800"
                "\t760\t800\t60\ttp:A:P\n",
                "u\t1000\t0\t750\t+\tchr2\t100000\t2000\t2750"
                "\t735\t750\t60\ttp:A:P\n",
            ]
        )
        result = evaluate_locus_evidence(
            records,
            unitig_ids=("u",),
            source_states={"u": "haplotig"},
            independent_read_support={},
            allowed_loci=("chr1", "chr2"),
            thresholds=self.thresholds,
            min_query_length=1,
            min_alignment_block_length=1,
            max_query_gap=1000,
            max_target_gap=1000,
        )
        decision = result.decisions[0]
        self.assertEqual(decision.status, "ambiguous")
        self.assertEqual(decision.reason, "competing_mT2T_loci")
        self.assertAlmostEqual(decision.next_best_margin or 0.0, 0.025)
        self.assertEqual(result.filtered_records, ())

    def test_unique_locus_filters_paf_without_assigning_a_group(self) -> None:
        records = parse_paf_lines(
            [
                "u\t1000\t0\t800\t+\tchr1\t100000\t1000\t1800"
                "\t760\t800\t60\tcg:Z:800M\ttp:A:P\n",
                "u\t1000\t0\t400\t+\tchr2\t100000\t2000\t2400"
                "\t320\t400\t60\ttp:A:P\n",
            ]
        )
        result = evaluate_locus_evidence(
            records,
            unitig_ids=("u",),
            source_states={"u": "haplotig"},
            independent_read_support={},
            allowed_loci=("chr1", "chr2"),
            thresholds=self.thresholds,
            min_query_length=1,
            min_alignment_block_length=1,
            max_query_gap=1000,
            max_target_gap=1000,
        )
        decision = result.decisions[0]
        self.assertEqual(decision.status, "locus_assigned_haplotype_unresolved")
        self.assertEqual(decision.assigned_locus, "chr1")
        self.assertIsNone(decision.assigned_group)
        self.assertEqual(
            tuple(record.target_name for record in result.filtered_records),
            ("chr1",),
        )

    def test_low_coverage_needs_independent_support_even_with_unique_locus(self) -> None:
        records = parse_paf_lines(
            [
                "low\t1000\t0\t800\t+\tchr1\t100000\t1000\t1800"
                "\t760\t800\t60\ttp:A:P\n",
            ]
        )
        without_support = evaluate_locus_evidence(
            records,
            unitig_ids=("low",),
            source_states={"low": "low_coverage"},
            independent_read_support={},
            allowed_loci=("chr1",),
            thresholds=self.thresholds,
            min_query_length=1,
            min_alignment_block_length=1,
            max_query_gap=1000,
            max_target_gap=1000,
        )
        with_support = evaluate_locus_evidence(
            records,
            unitig_ids=("low",),
            source_states={"low": "low_coverage"},
            independent_read_support={"low": 2},
            allowed_loci=("chr1",),
            thresholds=self.thresholds,
            min_query_length=1,
            min_alignment_block_length=1,
            max_query_gap=1000,
            max_target_gap=1000,
        )
        self.assertEqual(without_support.decisions[0].status, "unassigned")
        self.assertEqual(without_support.filtered_records, ())
        self.assertEqual(
            with_support.decisions[0].status,
            "locus_assigned_haplotype_unresolved",
        )
        self.assertEqual(
            with_support.decisions[0].rescue_class,
            "low_coverage_haplotig_candidate",
        )

    def test_secondary_and_missing_tp_are_audited_not_used(self) -> None:
        records = parse_paf_lines(
            [
                "u\t1000\t0\t800\t+\tchr1\t100000\t0\t800"
                "\t800\t800\t60\ttp:A:S\n",
                "u\t1000\t0\t800\t+\tchr2\t100000\t0\t800"
                "\t800\t800\t60\n",
            ]
        )
        result = evaluate_locus_evidence(
            records,
            unitig_ids=("u",),
            source_states={"u": "haplotig"},
            independent_read_support={},
            allowed_loci=("chr1", "chr2"),
            thresholds=self.thresholds,
            min_query_length=1,
            min_alignment_block_length=1,
            max_query_gap=1000,
            max_target_gap=1000,
        )
        self.assertEqual(
            tuple(row.reason for row in result.alignment_audit),
            ("secondary_alignment", "missing_tp_tag"),
        )
        self.assertEqual(result.decisions[0].status, "unassigned")

    def test_duplicate_alignment_is_audited_and_counted_once(self) -> None:
        records = parse_paf_lines(
            [
                "u\t1000\t0\t800\t+\tchr1\t100000\t0\t800"
                "\t800\t800\t60\ttp:A:P\n",
                "u\t1000\t0\t800\t+\tchr1\t100000\t0\t800"
                "\t800\t800\t60\tcg:Z:800M\ttp:A:P\n",
            ]
        )
        result = evaluate_locus_evidence(
            records,
            unitig_ids=("u",),
            source_states={"u": "haplotig"},
            independent_read_support={},
            allowed_loci=("chr1",),
            thresholds=self.thresholds,
            min_query_length=1,
            min_alignment_block_length=1,
            max_query_gap=1000,
            max_target_gap=1000,
        )
        self.assertEqual(
            tuple(row.reason for row in result.alignment_audit),
            ("chain_candidate", "duplicate_alignment"),
        )
        self.assertEqual(result.chains[0].record_count, 1)
        self.assertEqual(len(result.filtered_records), 1)

    def test_outputs_are_byte_stable_when_input_order_changes(self) -> None:
        lines = [
            "b\t1000\t0\t800\t+\tchrB\t100000\t0\t800"
            "\t800\t800\t60\ttp:A:P\n",
            "a\t1000\t0\t800\t+\tchrA\t100000\t0\t800"
            "\t800\t800\t60\ttp:A:P\n",
        ]
        outputs = []
        for ordered in (lines, list(reversed(lines))):
            records = parse_paf_lines(ordered, source="input.paf")
            targets = select_longest_targets(records, target_count=2)
            result = evaluate_locus_evidence(
                records,
                unitig_ids=("a", "b"),
                source_states={"a": "haplotig", "b": "haplotig"},
                independent_read_support={},
                allowed_loci=targets,
                thresholds=self.thresholds,
                min_query_length=1,
                min_alignment_block_length=1,
                max_query_gap=1000,
                max_target_gap=1000,
            )
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                paths = [root / f"out{index}.tsv" for index in range(5)]
                write_outputs(
                    result,
                    filtered_paf=str(paths[0]),
                    alignment_audit=str(paths[1]),
                    candidate_audit=str(paths[2]),
                    decision_audit=str(paths[3]),
                    routing_audit=str(paths[4]),
                )
                outputs.append(tuple(path.read_bytes() for path in paths))
        stable_output_indexes = (0, 2, 3, 4)
        self.assertEqual(
            tuple(outputs[0][index] for index in stable_output_indexes),
            tuple(outputs[1][index] for index in stable_output_indexes),
        )


if __name__ == "__main__":
    unittest.main()
