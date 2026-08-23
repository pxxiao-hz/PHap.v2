import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CHAIN_SCRIPT = ROOT / "utils" / "find_collinear_chains.py"
TABLE_SCRIPT = ROOT / "utils" / "allelic_table_generate_v2.py"


def paf_line(
    query,
    query_length,
    query_start,
    query_end,
    strand,
    target,
    target_start,
    target_end,
    identity=0.98,
    target_length=1_000_000,
    mapq=60,
):
    block_length = query_end - query_start
    matches = round(block_length * identity)
    return "\t".join(
        map(
            str,
            [
                query,
                query_length,
                query_start,
                query_end,
                strand,
                target,
                target_length,
                target_start,
                target_end,
                matches,
                block_length,
                mapq,
                "tp:A:P",
            ],
        )
    )


class CollinearChainTests(unittest.TestCase):
    def test_selects_all_local_blocks_on_best_target_and_rejects_ambiguous_reference(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "raw.paf"
            output = work / "chain.paf"
            qc = work / "chain.qc.tsv"
            records = [
                        paf_line("q1", 100_000, 0, 40_000, "+", "chr1", 0, 40_000),
                        paf_line("q1", 100_000, 40_000, 80_000, "+", "chr1", 40_000, 80_000),
                        # Query order increases but target order moves backwards.
                        paf_line("q1", 100_000, 80_000, 100_000, "+", "chr1", 10_000, 30_000),
                        paf_line("q1", 100_000, 0, 70_000, "+", "chr2", 0, 70_000, 0.97),
                        paf_line("q2", 100_000, 0, 80_000, "+", "chr1", 200_000, 280_000),
                        paf_line("q2", 100_000, 0, 80_000, "+", "chr2", 200_000, 280_000),
                        # Collinear but too sparse across its target span.
                        paf_line("q3", 200_000, 0, 20_000, "+", "chr3", 0, 20_000),
                        paf_line("q3", 200_000, 20_000, 40_000, "+", "chr3", 140_000, 160_000),
                        paf_line("q3", 200_000, 40_000, 60_000, "+", "chr3", 280_000, 300_000),
            ]
            records.extend(
                paf_line(
                    "q4",
                    5_000_000,
                    index * 250_000,
                    index * 250_000 + 50_000,
                    "+",
                    "chr4",
                    index * 250_000,
                    index * 250_000 + 50_000,
                    target_length=5_000_000,
                )
                for index in range(20)
            )
            records.extend(
                paf_line(
                    "q5",
                    10_000_000,
                    index * 200_000,
                    index * 200_000 + 25_000,
                    "+",
                    "chr5",
                    index * 200_000,
                    index * 200_000 + 25_000,
                    target_length=10_000_000,
                )
                for index in range(20)
            )
            paf.write_text("\n".join(records) + "\n")

            subprocess.run(
                [
                    sys.executable,
                    str(CHAIN_SCRIPT),
                    "--paf",
                    str(paf),
                    "--output",
                    str(output),
                    "--qc",
                    str(qc),
                    "--max-overlap",
                    "1000",
                    "--max-gap",
                    "200000",
                ],
                check=True,
            )

            selected = [line.split("\t") for line in output.read_text().splitlines()]
            q1_selected = [row for row in selected if row[0] == "q1"]
            self.assertEqual([(row[0], row[5]) for row in q1_selected], [("q1", "chr1")] * 3)
            self.assertEqual(
                [(int(row[2]), int(row[3])) for row in q1_selected],
                [(0, 40_000), (40_000, 80_000), (80_000, 100_000)],
            )

            with qc.open() as handle:
                rows = {row["query"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["q1"]["status"], "accepted")
            self.assertEqual(rows["q1"]["query_covered_bp"], "100000")
            self.assertEqual(rows["q2"]["reason"], "ambiguous_reference_chromosome")
            self.assertEqual(rows["q3"]["status"], "accepted")
            self.assertEqual(rows["q4"]["status"], "accepted")
            self.assertEqual(rows["q4"]["placement_mode"], "segmented")
            self.assertTrue(
                all("pv:Z:segmented" in row[12:] for row in selected if row[0] == "q4")
            )
            self.assertEqual(rows["q5"]["status"], "accepted")
            self.assertEqual(rows["q5"]["placement_mode"], "segmented")
            self.assertTrue(
                all(
                    "pv:Z:segmented" in row[12:]
                    for row in selected if row[0] == "q5"
                )
            )

    def test_accepts_mixed_strands_and_direction_changes_on_one_chromosome(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "raw.paf"
            output = work / "selected.paf"
            qc = work / "selection.qc.tsv"
            paf.write_text(
                "\n".join(
                    [
                        paf_line("mixed", 400_000, 0, 100_000, "+", "chr7", 0, 100_000),
                        paf_line("mixed", 400_000, 100_000, 200_000, "-", "chr7", 300_000, 400_000),
                        paf_line("mixed", 400_000, 200_000, 300_000, "+", "chr7", 100_000, 200_000),
                    ]
                )
                + "\n"
            )
            subprocess.run(
                [
                    sys.executable,
                    str(CHAIN_SCRIPT),
                    "--paf",
                    str(paf),
                    "--output",
                    str(output),
                    "--qc",
                    str(qc),
                ],
                check=True,
            )

            self.assertEqual(len(output.read_text().splitlines()), 3)
            with qc.open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["status"], "accepted")
            self.assertEqual(row["best_target"], "chr7")
            self.assertEqual(row["best_strand"], "mixed")
            self.assertEqual(row["strand_switches"], "2")

    def test_does_not_reject_haplotype_divergence_after_local_identity_filter(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "raw.paf"
            output = work / "selected.paf"
            qc = work / "selection.qc.tsv"
            paf.write_text(
                paf_line(
                    "divergent",
                    200_000,
                    0,
                    150_000,
                    "+",
                    "chr9",
                    0,
                    150_000,
                    identity=0.91,
                )
                + "\n"
            )
            subprocess.run(
                [
                    sys.executable,
                    str(CHAIN_SCRIPT),
                    "--paf",
                    str(paf),
                    "--output",
                    str(output),
                    "--qc",
                    str(qc),
                ],
                check=True,
            )

            with qc.open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["status"], "accepted")
            self.assertEqual(len(output.read_text().splitlines()), 1)

    def test_minus_strand_chain_uses_oriented_target_coordinates(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "raw.paf"
            output = work / "chain.paf"
            qc = work / "chain.qc.tsv"
            paf.write_text(
                "\n".join(
                    [
                        paf_line("minus", 100_000, 0, 40_000, "-", "chr1", 900_000, 940_000),
                        paf_line("minus", 100_000, 40_000, 80_000, "-", "chr1", 860_000, 900_000),
                    ]
                )
                + "\n"
            )
            subprocess.run(
                [
                    sys.executable,
                    str(CHAIN_SCRIPT),
                    "--paf",
                    str(paf),
                    "--output",
                    str(output),
                    "--qc",
                    str(qc),
                ],
                check=True,
            )
            self.assertEqual(len(output.read_text().splitlines()), 2)


class AllelicTableTests(unittest.TestCase):
    def write_table_inputs(self, work):
        paf = work / "chain.paf"
        contig_types = work / "contig_type.txt"
        records = [
                    paf_line("A", 100_000, 0, 100_000, "+", "chr1", 0, 100_000),
                    paf_line("C", 100_000, 0, 100_000, "+", "chr1", 0, 100_000),
                    paf_line("D", 100_000, 0, 100_000, "+", "chr1", 0, 100_000),
                    paf_line("B", 100_000, 0, 100_000, "+", "chr1", 100_000, 200_000),
                    # Equal-confidence alternatives exceed tetraploid capacity.
                    paf_line("E", 100_000, 0, 100_000, "+", "chr2", 0, 100_000),
                    paf_line("F", 100_000, 0, 100_000, "+", "chr2", 0, 100_000),
                    paf_line("G", 100_000, 0, 100_000, "+", "chr2", 0, 100_000),
                    paf_line("H", 100_000, 0, 100_000, "+", "chr2", 0, 100_000),
        ]
        records.extend(
            paf_line(
                "I",
                110_000,
                index * 10_000,
                (index + 1) * 10_000,
                "+",
                "chr3",
                index * 40_000,
                index * 40_000 + 10_000,
            )
            for index in range(11)
        )
        paf.write_text("\n".join(records) + "\n")
        contig_types.write_text(
            "contig_ID\tcontig_type\n"
            "A\thaplotig\nC\thaplotig\nD\tdiplotig\nB\thaplotig\n"
            "E\thaplotig\nF\thaplotig\nG\thaplotig\nH\tdiplotig\nI\thaplotig\n"
        )
        return paf, contig_types

    def run_table(self, work, paf, contig_types, gfa=None, extra_args=None):
        paths = {
            "output": work / "table.txt",
            "projections": work / "projections.tsv",
            "qc": work / "qc.tsv",
            "rejected": work / "rejected.tsv",
            "pairs": work / "pairs.tsv",
            "summary": work / "summary.json",
        }
        command = [
            sys.executable,
            str(TABLE_SCRIPT),
            "--paf",
            str(paf),
            "--contig-type",
            str(contig_types),
        ]
        for option, path in paths.items():
            command.extend(["--" + option, str(path)])
        if gfa is not None:
            command.extend(["--gfa", str(gfa)])
        if extra_args:
            command.extend(extra_args)
        subprocess.run(command, check=True)
        return paths

    def test_real_overlap_dosage_and_ambiguous_over_capacity(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf, contig_types = self.write_table_inputs(work)
            paths = self.run_table(work, paf, contig_types)

            rows = [line.split("\t") for line in paths["output"].read_text().splitlines()]
            self.assertEqual([row for row in rows if row[0] == "chr1"], [
                ["chr1", "0", "100000", "A", "C", "D"],
                ["chr1", "100000", "200000", "B"],
            ])
            self.assertEqual(len([row for row in rows if row[0] == "chr3"]), 11)
            self.assertFalse(any("A" in row[3:] and "B" in row[3:] for row in rows))

            with paths["qc"].open() as handle:
                qc_rows = list(csv.DictReader(handle, delimiter="\t"))
            chr2 = [row for row in qc_rows if row["target"] == "chr2"]
            self.assertEqual(chr2[0]["status"], "over_capacity_ambiguous")
            self.assertEqual(chr2[0]["selected"], "")

            summary = json.loads(paths["summary"].read_text())
            self.assertEqual(summary["status_counts"]["over_capacity_ambiguous"], 1)
            with paths["rejected"].open() as handle:
                rejected = list(csv.DictReader(handle, delimiter="\t"))
            self.assertFalse(any(row["unitig"] == "I" for row in rejected))
            self.assertEqual(summary["fragmented_unitigs_above_qc_limit"], 1)

    def test_fragmented_unitig_keeps_dominant_block_with_block_level_metrics(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "chain.paf"
            contig_types = work / "contig_type.txt"
            records = [
                paf_line(
                    "long",
                    2_000_000,
                    0,
                    1_000_000,
                    "+",
                    "chr1",
                    0,
                    1_000_000,
                    identity=0.999,
                    target_length=2_000_000,
                )
            ]
            records.extend(
                paf_line(
                    "long",
                    2_000_000,
                    1_100_000 + index * 20_000,
                    1_110_000 + index * 20_000,
                    "+",
                    "chr1",
                    1_100_000 + index * 40_000,
                    1_110_000 + index * 40_000,
                    identity=0.85,
                    target_length=2_000_000,
                )
                for index in range(11)
            )
            paf.write_text("\n".join(records) + "\n")
            contig_types.write_text(
                "contig_ID\tcontig_type\nlong\thaplotig\n"
            )
            paths = self.run_table(work, paf, contig_types)

            with paths["projections"].open() as handle:
                projections = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(projections), 1)
            projection = projections[0]
            self.assertEqual(projection["block_count"], "12")
            self.assertEqual(projection["projection_class"], "anchor")
            self.assertAlmostEqual(float(projection["identity"]), 0.999, places=3)
            self.assertGreater(float(projection["target_aligned_fraction"]), 0.90)

            with paths["rejected"].open() as handle:
                rejected = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                sum(row["reason"] == "block_identity_below_threshold" for row in rejected),
                11,
            )
            summary = json.loads(paths["summary"].read_text())
            self.assertEqual(summary["fragmented_unitigs_above_qc_limit"], 1)

    def test_sparse_syntenic_chain_projects_its_supported_span(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "chain.paf"
            contig_types = work / "contig_type.txt"
            records = [
                paf_line(
                    "long_haplotype",
                    5_000_000,
                    index * 250_000,
                    index * 250_000 + 50_000,
                    "+",
                    "chr4",
                    index * 250_000,
                    index * 250_000 + 50_000,
                    target_length=5_000_000,
                )
                + "\tpv:Z:sparse_syntenic"
                for index in range(20)
            ]
            paf.write_text("\n".join(records) + "\n")
            contig_types.write_text(
                "contig_ID\tcontig_type\nlong_haplotype\thaplotig\n"
            )
            paths = self.run_table(work, paf, contig_types)

            self.assertEqual(
                paths["output"].read_text().splitlines(),
                ["chr4\t0\t4800000\tlong_haplotype"],
            )
            with paths["projections"].open() as handle:
                projection = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(projection["placement_mode"], "sparse_syntenic")
            self.assertLess(float(projection["projection_coverage"]), 0.30)

    def test_fragmented_syntenic_chain_does_not_bridge_unsupported_gaps(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "chain.paf"
            contig_types = work / "contig_type.txt"
            records = [
                paf_line(
                    "fragmented_haplotype",
                    1_000_000,
                    query_start,
                    query_start + 100_000,
                    "+",
                    "chr5",
                    target_start,
                    target_start + 100_000,
                    target_length=1_000_000,
                )
                + "\tpv:Z:fragmented_syntenic"
                for query_start, target_start in [(0, 0), (500_000, 500_000)]
            ]
            paf.write_text("\n".join(records) + "\n")
            contig_types.write_text(
                "contig_ID\tcontig_type\nfragmented_haplotype\thaplotig\n"
            )
            paths = self.run_table(work, paf, contig_types)

            self.assertEqual(
                paths["output"].read_text().splitlines(),
                [
                    "chr5\t0\t100000\tfragmented_haplotype",
                    "chr5\t500000\t600000\tfragmented_haplotype",
                ],
            )
            with paths["projections"].open() as handle:
                projections = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(projections), 2)
            self.assertTrue(
                all(row["placement_mode"] == "fragmented_syntenic" for row in projections)
            )

    def test_segmented_projection_accepts_mixed_strands_and_nonmonotonic_blocks(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "selected.paf"
            contig_types = work / "contig_type.txt"
            records = [
                paf_line("rearranged", 1_000_000, 800_000, 900_000, "+", "chr7", 0, 100_000),
                paf_line("rearranged", 1_000_000, 400_000, 500_000, "+", "chr7", 100_000, 200_000),
                paf_line("rearranged", 1_000_000, 0, 100_000, "-", "chr7", 200_000, 300_000),
            ]
            paf.write_text(
                "\n".join(record + "\tpv:Z:segmented" for record in records) + "\n"
            )
            contig_types.write_text(
                "contig_ID\tcontig_type\nrearranged\thaplotig\n"
            )
            paths = self.run_table(work, paf, contig_types)

            self.assertEqual(
                paths["output"].read_text().splitlines(),
                ["chr7\t0\t300000\trearranged"],
            )
            with paths["projections"].open() as handle:
                projection = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(projection["placement_mode"], "segmented")
            self.assertEqual(projection["strand"], "mixed")
            self.assertLess(float(projection["collinearity"]), 0.80)

    def test_segmented_sparse_block_is_split_instead_of_rejected(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "selected.paf"
            contig_types = work / "contig_type.txt"
            records = [
                paf_line("split", 100_000, 0, 15_000, "+", "chr7", 0, 15_000),
                paf_line("split", 100_000, 15_000, 30_000, "+", "chr7", 35_000, 50_000),
            ]
            paf.write_text(
                "\n".join(record + "\tpv:Z:segmented" for record in records) + "\n"
            )
            contig_types.write_text("contig_ID\tcontig_type\nsplit\thaplotig\n")
            paths = self.run_table(work, paf, contig_types)

            self.assertEqual(
                paths["output"].read_text().splitlines(),
                ["chr7\t0\t15000\tsplit", "chr7\t35000\t50000\tsplit"],
            )
            with paths["projections"].open() as handle:
                projections = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(projections), 2)
            self.assertTrue(
                all(float(row["projection_coverage"]) == 1.0 for row in projections)
            )

    def test_long_segmented_path_envelope_excludes_lower_priority_fragments(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "selected.paf"
            contig_types = work / "contig_type.txt"
            records = [
                paf_line("long", 10_000_000, 0, 2_000_000, "+", "chr1", 0, 2_000_000),
                paf_line(
                    "long", 10_000_000, 4_000_000, 6_000_000, "+", "chr1", 4_000_000, 6_000_000
                ),
                paf_line(
                    "long", 10_000_000, 8_000_000, 10_000_000, "+", "chr1", 8_000_000, 10_000_000
                ),
                paf_line("fragment_a", 1_000_000, 0, 1_000_000, "+", "chr1", 2_200_000, 3_200_000),
                paf_line("fragment_b", 1_000_000, 0, 1_000_000, "+", "chr1", 2_200_000, 3_200_000),
            ]
            paf.write_text(
                "\n".join(
                    record + ("\tpv:Z:segmented" if record.startswith("long\t") else "")
                    for record in records
                )
                + "\n"
            )
            contig_types.write_text(
                "contig_ID\tcontig_type\n"
                "long\thaplotig\n"
                "fragment_a\tdiplotig\n"
                "fragment_b\tdiplotig\n"
            )
            paths = self.run_table(
                work,
                paf,
                contig_types,
                extra_args=[
                    "--max-projection-blocks",
                    "2",
                    "--over-capacity-policy",
                    "best",
                ],
            )

            with paths["projections"].open() as handle:
                projections = list(csv.DictReader(handle, delimiter="\t"))
            long_projections = [row for row in projections if row["unitig"] == "long"]
            self.assertEqual(
                sum(row["constraint_role"] == "path_envelope" for row in long_projections),
                1,
            )
            self.assertEqual(
                sum(row["constraint_role"] == "evidence_block" for row in long_projections),
                3,
            )

            with paths["pairs"].open() as handle:
                pairs = list(csv.DictReader(handle, delimiter="\t"))
            emitted_pairs = {
                frozenset((row["unitig1"], row["unitig2"])) for row in pairs
            }
            retained_fragments = {
                fragment
                for fragment in ("fragment_a", "fragment_b")
                if frozenset(("long", fragment)) in emitted_pairs
            }
            self.assertEqual(len(retained_fragments), 1)
            with paths["output"].with_name(
                "deferred_overcapacity_unitigs.tsv"
            ).open() as handle:
                deferred = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(deferred), 1)
            self.assertNotIn(deferred[0]["unitig"], retained_fragments)

            with paths["qc"].open() as handle:
                qc = list(csv.DictReader(handle, delimiter="\t"))
            crowded = [
                row
                for row in qc
                if row["start"] == "2200000" and row["end"] == "3200000"
            ]
            self.assertEqual(len(crowded), 1)
            self.assertIn("long", crowded[0]["selected"].split(","))

            summary = json.loads(paths["summary"].read_text())
            self.assertEqual(summary["long_path_constraints"]["path_envelopes"], 1)
            self.assertGreaterEqual(
                summary["long_path_constraints"]["protected_long_path_pairs"], 1
            )
            self.assertEqual(
                summary["deferred_over_capacity_unitigs"]["count"], 1
            )

    def test_recovers_high_coverage_pairs_from_ambiguous_over_capacity_segment(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "selected.paf"
            contig_types = work / "contig_type.txt"
            names = ("anchor_a", "anchor_b", "anchor_c", "fragment", "collapsed")
            paf.write_text(
                "\n".join(
                    paf_line(name, 1_000_000, 0, 1_000_000, "+", "chr1", 0, 1_000_000)
                    for name in names
                )
                + "\n"
            )
            contig_types.write_text(
                "contig_ID\tcontig_type\n"
                "anchor_a\thaplotig\n"
                "anchor_b\thaplotig\n"
                "anchor_c\thaplotig\n"
                "fragment\thaplotig\n"
                "collapsed\tdiplotig\n"
            )
            paths = self.run_table(work, paf, contig_types)

            with paths["qc"].open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["status"], "over_capacity_ambiguous")
            self.assertEqual(row["selected"], "")

            with paths["pairs"].open() as handle:
                pairs = {
                    frozenset((row["unitig1"], row["unitig2"]))
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertIn(frozenset(("fragment", "collapsed")), pairs)

            summary = json.loads(paths["summary"].read_text())
            self.assertGreaterEqual(
                summary["over_capacity_pair_constraints"]["recovered_pairs"], 1
            )

    def test_direct_gfa_link_is_not_emitted_as_an_allelic_pair(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf = work / "chain.paf"
            contig_types = work / "contig_type.txt"
            gfa = work / "graph.gfa"
            paf.write_text(
                paf_line("left", 100_000, 0, 100_000, "+", "chr1", 0, 100_000)
                + "\n"
                + paf_line("right", 100_000, 0, 100_000, "+", "chr1", 0, 100_000)
                + "\n"
            )
            contig_types.write_text(
                "contig_ID\tcontig_type\nleft\thaplotig\nright\thaplotig\n"
            )
            gfa.write_text(
                "S\tleft\t*\tLN:i:100000\n"
                "S\tright\t*\tLN:i:100000\n"
                "L\tleft\t+\tright\t+\t20000M\n"
            )
            paths = self.run_table(work, paf, contig_types, gfa)

            self.assertEqual(paths["output"].read_text(), "")
            with paths["qc"].open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["status"], "graph_link_conflict_omitted")
            self.assertEqual(row["graph_link_conflicts"], "left|right")

    def test_public_cli_writes_manifest_and_legacy_table_alias(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            paf, contig_types = self.write_table_inputs(work)
            p_utg = work / "p_utg.fa"
            mt2t = work / "mt2t.fa"
            p_utg.write_text(">placeholder\nA\n")
            mt2t.write_text(">chr1\nA\n")
            output = work / "run"

            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "PHap.py"),
                    "allelic_table",
                    "--p_utg",
                    str(p_utg),
                    "--mT2T",
                    str(mt2t),
                    "--contig_type",
                    str(contig_types),
                    "--paf",
                    str(paf),
                    "--output-dir",
                    str(output),
                ],
                check=True,
            )

            table = output / "corrected_allelic_table.txt"
            self.assertEqual(table.read_bytes(), (output / "allelic.ctg.table.v2").read_bytes())
            manifest = json.loads((output / "run_manifest.json").read_text())
            self.assertEqual(manifest["phap_version"], "2.0.0")
            self.assertEqual(Path(manifest["inputs"]["paf"]["path"]), paf.resolve())


if __name__ == "__main__":
    unittest.main()
