from __future__ import annotations

import unittest
from typing import Any

from phap_core.read_assignment import (
    alignment_is_usable,
    assign_reads,
    build_unitig_candidates,
    canonical_read_id,
    group_assigned_reads,
    paired_fastq_patterns,
    parse_unitig_dosages,
    stable_destination,
)


class ReadAssignmentTests(unittest.TestCase):
    def setUp(self) -> None:
        self.group_unitigs = {
            "chr01_h1": {"u_h1", "u_d12", "u_tri"},
            "chr01_h2": {"u_h2", "u_d12", "u_tri"},
            "chr01_h3": {"u_tri"},
        }
        self.dosage = {
            "u_h1": 1,
            "u_h2": 1,
            "u_d12": 2,
            "u_tri": 3,
        }
        self.candidates = build_unitig_candidates(
            self.group_unitigs,
            self.dosage,
            ploidy=3,
        )

    def test_collapsed_partition_is_disjoint_exhaustive_and_reproducible(self) -> None:
        reads = {f"read{i:03d}": {"u_d12"} for i in range(100)}
        first = assign_reads(reads, self.candidates, modality="hifi", seed=17)
        reordered = assign_reads(
            dict(reversed(list(reads.items()))),
            self.candidates,
            modality="hifi",
            seed=17,
        )
        self.assertEqual(first, reordered)

        grouped = group_assigned_reads(first, self.group_unitigs)
        h1 = set(grouped["chr01_h1"])
        h2 = set(grouped["chr01_h2"])
        self.assertFalse(h1 & h2)
        self.assertEqual(h1 | h2, set(reads))
        self.assertTrue(h1)
        self.assertTrue(h2)

    def test_multimapping_is_resolved_or_reported_explicitly(self) -> None:
        assignments = assign_reads(
            {
                "compatible": {"u_d12", "u_h1"},
                "conflicting": {"u_h1", "u_h2"},
                "unknown": {"u_missing"},
                "unaligned": set(),
            },
            self.candidates,
            modality="ont",
            seed=22,
        )
        by_read = {assignment.read_id: assignment for assignment in assignments}

        self.assertEqual(by_read["compatible"].status, "assigned")
        self.assertEqual(by_read["compatible"].destination_group, "chr01_h1")
        self.assertEqual(by_read["compatible"].reason, "compatible_multimap_resolved")
        self.assertEqual(by_read["conflicting"].status, "ambiguous")
        self.assertEqual(by_read["conflicting"].reason, "incompatible_multimap")
        self.assertEqual(by_read["unknown"].status, "unassigned")
        self.assertEqual(by_read["unaligned"].reason, "no_usable_alignment")

    def test_dosage_group_count_mismatch_is_not_guessed(self) -> None:
        candidates = build_unitig_candidates(
            {"chr01_h1": {"bad_collapsed"}},
            {"bad_collapsed": 2},
            ploidy=3,
        )
        self.assertEqual(candidates["bad_collapsed"].status, "ambiguous")
        assignment = assign_reads(
            {"read1": {"bad_collapsed"}},
            candidates,
            modality="hifi",
            seed=1,
        )[0]
        self.assertEqual(assignment.status, "ambiguous")
        self.assertIsNone(assignment.destination_group)

    def test_current_dosage_output_preserves_invalid_calls(self) -> None:
        dosages, labels = parse_unitig_dosages(
            [
                "contig_ID\taverage_depth\tcontig_type\tdosage\tstatus\n",
                "u_hap\t20\thaplotig\t1\tassigned\n",
                "u_di\t40\tdiplotig\t2\tassigned\n",
                "u_mixed\t30\tambiguous\t.\tmixed\n",
                "u_six\t120\tdosage_6\t6\tassigned\n",
            ]
        )
        self.assertEqual(dosages["u_hap"], 1)
        self.assertEqual(dosages["u_di"], 2)
        self.assertIsNone(dosages["u_mixed"])
        self.assertEqual(dosages["u_six"], 6)
        self.assertEqual(labels["u_mixed"], "ambiguous")

        candidates = build_unitig_candidates(
            {"h1": {"u_mixed"}},
            dosages,
            ploidy=6,
        )
        self.assertEqual(candidates["u_mixed"].status, "unassigned")
        self.assertEqual(candidates["u_mixed"].reason, "invalid_dosage_call")

    def test_legacy_dosage_labels_remain_supported(self) -> None:
        dosages, _ = parse_unitig_dosages(
            [
                "legacy_hap\t20\thaplotig\n",
                "legacy_tri\t60\ttriplotig\n",
            ]
        )
        self.assertEqual(dosages, {"legacy_hap": 1, "legacy_tri": 3})

    def test_hash_uses_modality_and_not_python_hash_state(self) -> None:
        self.assertEqual(
            stable_destination(
                read_id="pair001",
                modality="hic",
                unitigs=["u_tri"],
                candidate_groups=["chr01_h1", "chr01_h2", "chr01_h3"],
                seed=99,
            ),
            stable_destination(
                read_id="pair001",
                modality="hic",
                unitigs=["u_tri"],
                candidate_groups=["chr01_h1", "chr01_h2", "chr01_h3"],
                seed=99,
            ),
        )
        observed = {
            stable_destination(
                read_id=f"pair{i}",
                modality="hic",
                unitigs=["u_tri"],
                candidate_groups=["chr01_h1", "chr01_h2", "chr01_h3"],
                seed=99,
            )
            for i in range(60)
        }
        self.assertEqual(observed, {"chr01_h1", "chr01_h2", "chr01_h3"})

    def test_ploidy_is_not_hard_coded_to_four(self) -> None:
        groups = {f"h{i}": {"u6"} for i in range(1, 7)}
        candidates = build_unitig_candidates(groups, {"u6": 6}, ploidy=6)
        self.assertEqual(candidates["u6"].status, "assignable")
        self.assertEqual(len(candidates["u6"].groups), 6)

    def test_hic_mate_suffixes_share_one_canonical_entity(self) -> None:
        self.assertEqual(canonical_read_id("pair001/1", paired=True), "pair001")
        self.assertEqual(canonical_read_id("pair001/2", paired=True), "pair001")
        self.assertEqual(canonical_read_id("longread/1", paired=False), "longread/1")
        self.assertEqual(
            paired_fastq_patterns("pair001/1", mate=1),
            ("pair001", "pair001/1"),
        )
        self.assertEqual(
            paired_fastq_patterns("pair001/2", mate=2),
            ("pair001", "pair001/2"),
        )

    def test_bam_filter_policy_keeps_nonproper_but_filters_bad_records(self) -> None:
        baseline: dict[str, Any] = {
            "is_unmapped": False,
            "is_secondary": False,
            "is_supplementary": False,
            "is_duplicate": False,
            "is_qcfail": False,
            "is_proper_pair": False,
            "mapping_quality": 20,
            "min_mapq": 10,
        }
        self.assertTrue(alignment_is_usable(**baseline))
        for flag in (
            "is_unmapped",
            "is_secondary",
            "is_supplementary",
            "is_duplicate",
            "is_qcfail",
        ):
            values = dict(baseline)
            values[flag] = True
            self.assertFalse(alignment_is_usable(**values), flag)
        low_mapq = dict(baseline)
        low_mapq["mapping_quality"] = 9
        self.assertFalse(alignment_is_usable(**low_mapq))


if __name__ == "__main__":
    unittest.main()
