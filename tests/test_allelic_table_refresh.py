from __future__ import annotations

import unittest

from phap_core.allelic_table import (
    AllelicTableRow,
    filter_allelic_rows,
    find_missing_bridge_unitigs,
)


class AllelicTableRefreshTests(unittest.TestCase):
    def test_low_coverage_is_excluded_with_audit_not_coerced_to_haplotig(self) -> None:
        rows, audit = filter_allelic_rows(
            [AllelicTableRow("chr1", 0, 100, ("hap", "low"))],
            {"hap": 1, "low": None},
            {"hap": "haplotig", "low": "low_coverage"},
            {"hap": 1000, "low": 900},
            ploidy=4,
        )
        self.assertEqual(rows, (AllelicTableRow("chr1", 0, 100, ("hap",)),))
        self.assertEqual(audit[0].unitig_id, "low")
        self.assertEqual(audit[0].source_state, "low_coverage")
        self.assertEqual(audit[0].reason, "invalid_or_missing_dosage")

    def test_ploidy_six_capacity_is_not_hard_coded_to_four(self) -> None:
        rows, audit = filter_allelic_rows(
            [AllelicTableRow("chr1", 0, 100, ("u3", "u2", "u1"))],
            {"u3": 3, "u2": 2, "u1": 1},
            {"u3": "dosage_3", "u2": "dosage_2", "u1": "haplotig"},
            {"u3": 3000, "u2": 2000, "u1": 1000},
            ploidy=6,
        )
        self.assertEqual(
            rows,
            (AllelicTableRow("chr1", 0, 100, ("u1", "u2", "u3")),),
        )
        self.assertEqual(audit, ())

    def test_bin_with_only_invalid_dosage_records_is_not_written_empty(self) -> None:
        rows, audit = filter_allelic_rows(
            [AllelicTableRow("chr1", 0, 100, ("low",))],
            {"low": None},
            {"low": "low_coverage"},
            {"low": 900},
            ploidy=4,
        )
        self.assertEqual(rows, (None,))
        self.assertEqual(len(audit), 1)
        self.assertEqual(audit[0].unitig_id, "low")
        self.assertEqual(audit[0].reason, "invalid_or_missing_dosage")

    def test_unresolvable_over_capacity_bin_is_not_partially_guessed(self) -> None:
        rows, audit = filter_allelic_rows(
            [AllelicTableRow("chr1", 0, 100, ("u3", "u2"))],
            {"u3": 3, "u2": 2},
            {"u3": "triplotig", "u2": "diplotig"},
            {"u3": 3000, "u2": 2000},
            ploidy=4,
        )
        self.assertEqual(rows, (None,))
        self.assertTrue(
            all(row.reason == "ambiguous_bin_exceeds_ploidy" for row in audit)
        )

    def test_refresh_window_does_not_cross_chromosome_boundaries(self) -> None:
        rows = [
            AllelicTableRow("chr1", 0, 100, ("cross",)),
            AllelicTableRow("chr2", 0, 100, ("local",)),
            AllelicTableRow("chr2", 100, 200, ("cross",)),
        ]
        missing = find_missing_bridge_unitigs(
            rows,
            1,
            search_range=1,
        )
        self.assertEqual(missing, ())

    def test_refresh_fills_a_same_chromosome_internal_gap(self) -> None:
        rows = [
            AllelicTableRow("chr1", 0, 100, ("bridge",)),
            AllelicTableRow("chr1", 100, 200, ("local",)),
            AllelicTableRow("chr1", 200, 300, ("bridge",)),
        ]
        missing = find_missing_bridge_unitigs(
            rows,
            1,
            search_range=1,
        )
        self.assertEqual(missing, ("bridge",))


if __name__ == "__main__":
    unittest.main()
