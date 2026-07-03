from __future__ import annotations

import unittest

from phap_core.paf import (
    PafFormatError,
    parse_paf_line,
    parse_paf_lines,
    select_primary_records,
)


class PafParsingTests(unittest.TestCase):
    def test_optional_tags_are_looked_up_by_name_not_position(self) -> None:
        first = parse_paf_line(
            "u1\t1000\t0\t800\t+\tchr1\t10000\t10\t810\t780\t800\t60"
            "\ttp:A:P\tcg:Z:800M\tNM:i:20"
        )
        second = parse_paf_line(
            "u2\t1000\t0\t800\t+\tchr1\t10000\t10\t810\t780\t800\t60"
            "\tcg:Z:800M\tNM:i:20\ttp:A:P"
        )
        self.assertTrue(first.is_primary)
        self.assertTrue(second.is_primary)
        self.assertEqual(first.tag("tp"), first.tags[0])
        self.assertEqual(second.tag("tp"), second.tags[2])

    def test_optional_tag_value_can_contain_colons(self) -> None:
        record = parse_paf_line(
            "u1\t100\t0\t8\t+\tchr1\t1000\t0\t8\t7\t8\t60"
            "\tcs:Z::5*ac:2\ttp:A:P"
        )
        cs_tag = record.tag("cs")
        self.assertIsNotNone(cs_tag)
        self.assertEqual(cs_tag.value if cs_tag is not None else None, ":5*ac:2")
        self.assertEqual(
            record.to_line(),
            "u1\t100\t0\t8\t+\tchr1\t1000\t0\t8\t7\t8\t60"
            "\tcs:Z::5*ac:2\ttp:A:P",
        )

    def test_primary_policy_audits_secondary_and_missing_tp(self) -> None:
        records = parse_paf_lines(
            [
                "p\t100\t0\t80\t+\tchr1\t1000\t0\t80\t80\t80\t60\ttp:A:P\n",
                "i1\t100\t0\t80\t-\tchr1\t1000\t0\t80\t80\t80\t60\ttp:A:I\n",
                "s\t100\t0\t80\t+\tchr1\t1000\t0\t80\t80\t80\t60\ttp:A:S\n",
                "i2\t100\t0\t80\t-\tchr1\t1000\t0\t80\t80\t80\t60\ttp:A:i\n",
                "m\t100\t0\t80\t+\tchr1\t1000\t0\t80\t80\t80\t60\n",
            ]
        )
        accepted, audit = select_primary_records(records, require_tp=True)
        self.assertEqual(
            tuple(record.query_name for record in accepted),
            ("i1", "p"),
        )
        self.assertEqual(
            {row.query_name: row.reason for row in audit},
            {
                "p": "primary_alignment",
                "i1": "primary_inversion",
                "s": "secondary_alignment",
                "i2": "secondary_inversion",
                "m": "missing_tp_tag",
            },
        )

    def test_reverse_record_keeps_forward_query_coordinates(self) -> None:
        record = parse_paf_line(
            "rev\t1000\t400\t900\t-\tchr3\t20000\t9100\t9600"
            "\t475\t500\t60\ttp:A:P"
        )
        self.assertEqual((record.query_start, record.query_end), (400, 900))
        self.assertEqual((record.target_start, record.target_end), (9100, 9600))
        self.assertEqual(record.strand, "-")

    def test_malformed_rows_report_source_and_line(self) -> None:
        malformed = (
            "u\t1000\t0\t100\t+\tchr1\t1000\t0\t100\t101\t100\t60"
        )
        with self.assertRaisesRegex(
            PafFormatError,
            r"fixture\.paf:7: matching bases",
        ):
            parse_paf_line(malformed, source="fixture.paf", line_number=7)

    def test_duplicate_tag_is_rejected(self) -> None:
        line = (
            "u\t1000\t0\t100\t+\tchr1\t1000\t0\t100\t100\t100\t60"
            "\ttp:A:P\ttp:A:S"
        )
        with self.assertRaisesRegex(PafFormatError, "duplicate optional tag"):
            parse_paf_line(line)

    def test_invalid_half_open_interval_is_rejected(self) -> None:
        line = "u\t1000\t0\t1001\t+\tchr1\t1000\t0\t100\t100\t100\t60"
        with self.assertRaisesRegex(PafFormatError, "half-open query interval"):
            parse_paf_line(line)


if __name__ == "__main__":
    unittest.main()
