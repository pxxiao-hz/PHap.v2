from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from typing import Optional

from phap_core.fasta import FastaRecord
from phap_core.mt2t_overlap import (
    Mt2tOverlapParameters,
    assemble_mt2t,
    reverse_complement,
)
from phap_core.paf import parse_paf_lines
from utils.mt2t_overlap_workflow import run_mt2t_overlap_workflow


def paf(
    query: str,
    query_length: int,
    query_start: int,
    query_end: int,
    strand: str,
    target: str,
    target_length: int,
    target_start: int,
    target_end: int,
    matching_bases: Optional[int] = None,
) -> str:
    block = query_end - query_start
    matching = block if matching_bases is None else matching_bases
    return (
        f"{query}\t{query_length}\t{query_start}\t{query_end}\t{strand}\t"
        f"{target}\t{target_length}\t{target_start}\t{target_end}\t"
        f"{matching}\t{block}\t60\ttp:A:P\n"
    )


class Mt2tOverlapTests(unittest.TestCase):
    def setUp(self) -> None:
        self.parameters = Mt2tOverlapParameters(
            min_contig_length=1,
            min_alignment_block_length=1,
            min_identity=0.8,
            min_shorter_coverage=0.1,
            min_longer_coverage=0.05,
            internal_margin_ratio=0.1,
            max_query_gap=100,
            max_target_gap=100,
        )

    def records(self, sequences: dict[str, str]) -> tuple[FastaRecord, ...]:
        return tuple(
            FastaRecord(identifier, identifier, sequence)
            for identifier, sequence in sequences.items()
        )

    def test_forward_suffix_prefix_join(self) -> None:
        fasta = self.records({"A": "AAAACCCC", "B": "CCCCGGGG"})
        records = parse_paf_lines([paf("B", 8, 0, 4, "+", "A", 8, 4, 8)])
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(len(result.paths), 1)
        self.assertEqual(result.outputs[0].sequence, "AAAACCCCGGGG")
        self.assertEqual(result.joins[0].right_trim_bases, 4)
        self.assertEqual(
            {row.source_id: row.status for row in result.routing},
            {"A": "merged", "B": "merged"},
        )

    def test_reverse_join_uses_both_original_interval_endpoints(self) -> None:
        fasta = self.records({"A": "AAAACCCC", "B": "TTTTGGGG"})
        records = parse_paf_lines([paf("B", 8, 4, 8, "-", "A", 8, 4, 8)])
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(result.outputs[0].sequence, "AAAACCCCAAAA")
        join = result.joins[0]
        self.assertEqual(join.right_orientation, "-")
        self.assertEqual((join.right_start, join.right_end), (0, 4))
        self.assertEqual(join.right_trim_bases, 4)

    def test_nonzero_terminal_overhang_is_preserved_by_rejecting_join(self) -> None:
        fasta = self.records({"A": "A" * 10, "B": "C" * 10})
        records = parse_paf_lines([paf("B", 10, 2, 6, "+", "A", 10, 6, 10)])
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(result.paths, ())
        self.assertEqual(
            {row.identifier: row.sequence for row in result.outputs},
            {"A": "A" * 10, "B": "C" * 10},
        )

    def test_unaligned_gap_inside_overlap_span_is_not_trimmed(self) -> None:
        fasta = self.records({"A": "AAAACCCC", "B": "CCTTCCGGGG"})
        records = parse_paf_lines(
            [
                paf("B", 10, 0, 2, "+", "A", 8, 4, 6),
                paf("B", 10, 4, 6, "+", "A", 8, 6, 8),
            ]
        )
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(result.paths, ())
        self.assertEqual(
            {row.identifier: row.sequence for row in result.outputs},
            {"A": "AAAACCCC", "B": "CCTTCCGGGG"},
        )
        chain = result.chain_audit[0]
        self.assertEqual(chain.query_union_bases, 4)
        self.assertEqual(chain.query_span_bases, 6)

    def test_canonical_reverse_traversal_transforms_both_intervals(self) -> None:
        fasta = self.records({"Z": "TTTTACGT", "A": "ACGTGGGA"})
        records = parse_paf_lines([paf("A", 8, 0, 4, "+", "Z", 8, 4, 8)])
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(
            result.outputs[0].sequence,
            reverse_complement("TTTTACGTGGGA"),
        )
        join = result.joins[0]
        self.assertEqual(join.left_id, "A")
        self.assertEqual(join.left_orientation, "-")
        self.assertEqual(join.right_orientation, "-")
        self.assertEqual((join.right_start, join.right_end), (0, 4))

    def test_three_contig_path_consumes_each_contig_once(self) -> None:
        fasta = self.records(
            {"A": "AAAACCCC", "B": "CCCCGGGG", "C": "GGGGTTTT"}
        )
        records = parse_paf_lines(
            [
                paf("B", 8, 0, 4, "+", "A", 8, 4, 8),
                paf("C", 8, 0, 4, "+", "B", 8, 4, 8),
            ]
        )
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(result.outputs[0].sequence, "AAAACCCCGGGGTTTT")
        self.assertEqual(result.outputs[0].source_ids, ("A", "B", "C"))
        self.assertEqual(len({row.source_id for row in result.routing}), 3)

    def test_branch_component_is_preserved_and_audited(self) -> None:
        fasta = self.records(
            {"A": "AAAACCCC", "B": "CCCCGGGG", "C": "CCCCTTTT"}
        )
        records = parse_paf_lines(
            [
                paf("B", 8, 0, 4, "+", "A", 8, 4, 8),
                paf("C", 8, 0, 4, "+", "A", 8, 4, 8),
            ]
        )
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(result.paths, ())
        self.assertEqual({row.identifier for row in result.outputs}, {"A", "B", "C"})
        self.assertEqual(
            {row.reason for row in result.routing},
            {"branch_conflict"},
        )
        self.assertEqual(
            {row.reason for row in result.edge_audit},
            {"branch_conflict"},
        )

    def test_cycle_component_is_preserved(self) -> None:
        fasta = self.records(
            {"A": "AAAACCCC", "B": "CCCCGGGG", "C": "GGGGAAAA"}
        )
        records = parse_paf_lines(
            [
                paf("B", 8, 0, 4, "+", "A", 8, 4, 8),
                paf("C", 8, 0, 4, "+", "B", 8, 4, 8),
                paf("A", 8, 0, 4, "+", "C", 8, 4, 8),
            ]
        )
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(result.paths, ())
        self.assertEqual({row.reason for row in result.routing}, {"cycle_conflict"})

    def test_conflicting_pair_orientations_are_not_guessed(self) -> None:
        fasta = self.records({"A": "AAAACCCC", "B": "CCCCGGGG"})
        records = parse_paf_lines(
            [
                paf("B", 8, 0, 4, "+", "A", 8, 4, 8),
                paf("A", 8, 4, 8, "-", "B", 8, 4, 8),
            ]
        )
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(result.paths, ())
        self.assertEqual(
            {row.reason for row in result.routing},
            {"orientation_conflict"},
        )

    def test_unique_internal_containment_routes_to_host(self) -> None:
        fasta = self.records({"host": "AAAACCCCGGGG", "child": "CCCC"})
        records = parse_paf_lines(
            [paf("child", 4, 0, 4, "+", "host", 12, 4, 8)]
        )
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(tuple(row.identifier for row in result.outputs), ("host",))
        routing = {row.source_id: row for row in result.routing}
        self.assertEqual(routing["child"].destination_id, "host")
        self.assertEqual(routing["child"].status, "contained")
        self.assertEqual(result.contained_ids, ("child",))

    def test_partial_containment_preserves_unique_child_ends(self) -> None:
        fasta = self.records({"host": "A" * 20, "child": "C" * 10})
        records = parse_paf_lines(
            [paf("child", 10, 1, 9, "+", "host", 20, 6, 14)]
        )
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(result.contained_ids, ())
        self.assertEqual({row.identifier for row in result.outputs}, {"host", "child"})

    def test_low_identity_terminal_overlap_is_rejected(self) -> None:
        fasta = self.records({"A": "AAAACCCC", "B": "CCCCGGGG"})
        records = parse_paf_lines(
            [paf("B", 8, 0, 4, "+", "A", 8, 4, 8, matching_bases=1)]
        )
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(result.paths, ())
        self.assertEqual(result.chain_audit[0].reason, "identity_below_threshold")

    def test_short_long_coverage_is_independent_of_paf_direction(self) -> None:
        fasta = self.records({"A": "A" * 100, "B": "A" * 60})
        forward = parse_paf_lines(
            [paf("B", 60, 0, 8, "+", "A", 100, 92, 100)]
        )
        reciprocal = parse_paf_lines(
            [paf("A", 100, 92, 100, "+", "B", 60, 0, 8)]
        )
        forward_result = assemble_mt2t(fasta, forward, self.parameters)
        reciprocal_result = assemble_mt2t(fasta, reciprocal, self.parameters)
        self.assertEqual(len(forward_result.paths), 1)
        self.assertEqual(len(reciprocal_result.paths), 1)
        self.assertEqual(
            forward_result.outputs[0].sequence,
            reciprocal_result.outputs[0].sequence,
        )

    def test_competing_containment_hosts_preserve_the_child(self) -> None:
        fasta = self.records(
            {
                "host1": "AAAACCCCGGGG",
                "host2": "TTTTCCCCAAAA",
                "child": "CCCC",
            }
        )
        records = parse_paf_lines(
            [
                paf("child", 4, 0, 4, "+", "host1", 12, 4, 8),
                paf("child", 4, 0, 4, "+", "host2", 12, 4, 8),
            ]
        )
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(result.contained_ids, ())
        self.assertEqual(
            {row.identifier for row in result.outputs},
            {"child", "host1", "host2"},
        )
        child = next(row for row in result.routing if row.source_id == "child")
        self.assertEqual(child.status, "preserved")
        self.assertEqual(child.reason, "competing_containment_hosts")

    def test_secondary_alignment_is_excluded_by_named_tp_tag(self) -> None:
        fasta = self.records({"A": "AAAACCCC", "B": "CCCCGGGG"})
        line = paf("B", 8, 0, 4, "+", "A", 8, 4, 8).rstrip("\n")
        records = parse_paf_lines([line.replace("tp:A:P", "cg:Z:4M\ttp:A:S") + "\n"])
        result = assemble_mt2t(fasta, records, self.parameters)
        self.assertEqual(result.paths, ())
        self.assertEqual(result.alignment_audit[0].initial_reason, "secondary_alignment")

    def test_paf_unknown_id_and_length_mismatch_fail_fast(self) -> None:
        fasta = self.records({"A": "AAAA", "B": "CCCC"})
        unknown = parse_paf_lines([paf("missing", 4, 0, 2, "+", "A", 4, 2, 4)])
        mismatch = parse_paf_lines([paf("B", 5, 0, 2, "+", "A", 4, 2, 4)])
        with self.assertRaisesRegex(ValueError, "absent from FASTA"):
            assemble_mt2t(fasta, unknown, self.parameters)
        with self.assertRaisesRegex(ValueError, "lengths disagree"):
            assemble_mt2t(fasta, mismatch, self.parameters)

    def test_alignment_audit_keys_final_reason_by_source_and_line(self) -> None:
        fasta = self.records(
            {
                "A": "AAAACCCC",
                "B": "CCCCGGGG",
                "C": "TTTTGGGG",
                "D": "GGGGAAAA",
            }
        )
        first = parse_paf_lines(
            [paf("B", 8, 0, 4, "+", "A", 8, 4, 8, matching_bases=1)],
            source="low.paf",
        )
        second = parse_paf_lines(
            [paf("D", 8, 0, 4, "+", "C", 8, 4, 8)],
            source="high.paf",
        )
        result = assemble_mt2t(fasta, first + second, self.parameters)
        reasons = {
            (row.source, row.line_number): row.final_reason
            for row in result.alignment_audit
        }
        self.assertEqual(reasons[("low.paf", 1)], "identity_below_threshold")
        self.assertEqual(reasons[("high.paf", 1)], "accepted_unambiguous_path")

    def test_overlapping_records_use_interval_union_coverage(self) -> None:
        fasta = self.records({"A": "A" * 20, "B": "A" * 10})
        records = parse_paf_lines(
            [
                paf("B", 10, 0, 6, "+", "A", 20, 10, 16),
                paf("B", 10, 4, 10, "+", "A", 20, 14, 20),
            ]
        )
        result = assemble_mt2t(fasta, records, self.parameters)
        chain = result.chain_audit[0]
        self.assertEqual(chain.query_coverage, 1.0)
        self.assertEqual(chain.target_coverage, 0.5)

    def test_short_contig_is_preserved_in_output_and_routing(self) -> None:
        parameters = Mt2tOverlapParameters(min_contig_length=10)
        fasta = self.records({"short": "ACGT"})
        result = assemble_mt2t(fasta, (), parameters)
        self.assertEqual(result.outputs[0].identifier, "short")
        self.assertEqual(result.routing[0].reason, "below_min_contig_length_preserved")

    def test_iupac_reverse_complement(self) -> None:
        self.assertEqual(reverse_complement("ACGTRYMKBDHVN"), "NBDHVMKRYACGT")

    def test_workflow_writes_complete_routing_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta_path = root / "input.fa"
            fasta_path.write_text(
                ">A source header\nAAAACCCC\n>B\nCCCCGGGG\n",
                encoding="utf-8",
            )
            paf_path = root / "input.paf"
            paf_path.write_text(
                paf("B", 8, 0, 4, "+", "A", 8, 4, 8),
                encoding="utf-8",
            )
            output = root / "out"
            run_mt2t_overlap_workflow(
                fasta_path=str(fasta_path),
                paf_path=str(paf_path),
                output_directory=str(output),
                parameters=self.parameters,
                tool_versions={"minimap2": "test"},
            )
            assembled = (output / "mT2T.fa").read_text(encoding="utf-8")
            routing = (output / "mt2t_sequence_routing.tsv").read_text(
                encoding="utf-8"
            )
            manifest = (output / "mt2t_manifest.tsv").read_text(encoding="utf-8")
        self.assertIn("AAAACCCCGGGG", assembled)
        self.assertIn("A\tPHap_mT2T_path_000001\tmerged", routing)
        self.assertIn("B\tPHap_mT2T_path_000001\tmerged", routing)
        self.assertIn("tool_minimap2\ttest", manifest)


if __name__ == "__main__":
    unittest.main()
