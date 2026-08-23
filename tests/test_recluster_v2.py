import csv
import json
import pickle
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RECLUSTER = ROOT / "utils" / "chr_uncluster_recluster.py"
SPLITTER = ROOT / "utils" / "split_clm_by_groups_v2.py"


class ReclusterTests(unittest.TestCase):
    def make_inputs(self, work):
        names = ["A", "B", "C", "D", "X", "Y", "R", "P", "N", "T"]
        fasta = work / "chr1.putg.fa"
        fasta.write_text("".join(f">{name}\n" + "GATC" * 25 + "\n" for name in names))
        types = {name: "haplotig" for name in names}
        types["Y"] = "diplotig"
        types["T"] = "tetraplotig"
        type_file = work / "types.tsv"
        type_file.write_text(
            "contig_ID\taverage_depth\tcontig_type\n"
            + "".join(f"{name}\t20\t{types[name]}\n" for name in names)
        )
        seeds = work / "groups.txt"
        seeds.write_text(
            "group1\t1\tA\n"
            "group2\t1\tB\n"
            "group3\t1\tC\n"
            "group4\t1\tD\n"
        )
        links = {
            tuple(sorted(("X", "B"))): 100,
            tuple(sorted(("Y", "A"))): 80,
            tuple(sorted(("Y", "C"))): 80,
            tuple(sorted(("Y", "B"))): 1,
            tuple(sorted(("Y", "D"))): 1,
            tuple(sorted(("R", "X"))): 100,
            tuple(sorted(("P", "A"))): 50,
            tuple(sorted(("P", "B"))): 50,
            tuple(sorted(("P", "X"))): 50,
            tuple(sorted(("P", "R"))): 50,
            tuple(sorted(("P", "Y"))): 100,
        }
        link_file = work / "links.pkl"
        with link_file.open("wb") as handle:
            pickle.dump(links, handle)
        return fasta, type_file, seeds, link_file

    def run_recluster(self, work, extra=None, check=True):
        fasta, types, seeds, links = self.make_inputs(work)
        output = work / "output"
        command = [
            sys.executable,
            str(RECLUSTER),
            "--fasta",
            str(fasta),
            "--contig-type",
            str(types),
            "--full-links",
            str(links),
            "--clusters-file",
            str(seeds),
            "--output-dir",
            str(output),
            "--min-group-margin",
            "0.2",
        ]
        if extra:
            command.extend(extra)
        result = subprocess.run(command, check=check, capture_output=True, text=True)
        return output, result

    def test_exact_dosage_propagation_and_audited_deferral(self):
        with tempfile.TemporaryDirectory() as temporary:
            output, result = self.run_recluster(Path(temporary))
            self.assertIn("assigned 8/10", result.stdout)
            with (output / "recluster_assignments.tsv").open() as handle:
                rows = {row["unitig"]: row for row in csv.DictReader(handle, delimiter="\t")}

            self.assertEqual(rows["X"]["groups"], "2")
            self.assertEqual(rows["R"]["groups"], "2")
            self.assertEqual(rows["R"]["assignment_round"], "2")
            self.assertEqual(rows["Y"]["groups"], "1,3")
            self.assertEqual(rows["Y"]["group_count"], rows["Y"]["dosage"])
            self.assertEqual(rows["Y"]["decision_adjusted_assigned_hic_links"], "80.000000")
            self.assertEqual(rows["T"]["groups"], "1,2,3,4")
            self.assertEqual(rows["T"]["assignment_basis"], "dosage_all_groups")
            self.assertEqual(rows["P"]["status"], "unassigned")
            self.assertEqual(rows["P"]["assignment_basis"], "low_group_margin")
            self.assertEqual(rows["N"]["assignment_basis"], "no_assigned_hic")

            summary = json.loads((output / "recluster_summary.json").read_text())
            self.assertEqual(summary["parameters"]["hic_link_normalization"], "dosage")
            self.assertEqual(summary["assigned"]["unitigs"], 8)
            self.assertEqual(summary["unassigned"]["unitigs"], 2)
            self.assertEqual(summary["round_assignment_counts"], [2, 1])
            self.assertEqual(summary["validation"]["violations"], 0)
            self.assertEqual(
                set((output / "unassigned_unitigs.txt").read_text().split()),
                {"P", "N"},
            )
            self.assertFalse((output / "full.links.txt").exists())

    def test_allelic_block_jointly_assigns_weak_member_by_constraint(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            names = ["A", "B", "C", "D", "X", "Y"]
            fasta = work / "chr1.putg.fa"
            fasta.write_text(
                "".join(f">{name}\n" + "GATC" * 25 + "\n" for name in names)
            )
            types = work / "types.tsv"
            types.write_text(
                "contig_ID\taverage_depth\tcontig_type\n"
                + "".join(f"{name}\t20\thaplotig\n" for name in names)
            )
            seeds = work / "groups.txt"
            seeds.write_text(
                "group1\t1\tA\ngroup2\t1\tB\n"
                "group3\t1\tC\ngroup4\t1\tD\n"
            )
            table = work / "table.tsv"
            table.write_text("chr1\t0\t100\tA\tD\tX\tY\n")
            links = work / "links.pkl"
            with links.open("wb") as handle:
                pickle.dump({("B", "X"): 100}, handle)
            output = work / "output"
            subprocess.run(
                [
                    sys.executable, str(RECLUSTER),
                    "--fasta", str(fasta), "--contig-type", str(types),
                    "--full-links", str(links), "--clusters-file", str(seeds),
                    "--allelic-table", str(table), "--output-dir", str(output),
                ],
                check=True, capture_output=True, text=True,
            )
            with (output / "recluster_assignments.tsv").open() as handle:
                rows = {
                    row["unitig"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(rows["X"]["groups"], "2")
            self.assertEqual(rows["Y"]["groups"], "3")
            self.assertEqual(rows["Y"]["assignment_basis"], "allelic_block_joint")
            self.assertEqual(
                rows["Y"]["allelic_block_acceptance_basis"], "constraint_forced"
            )
            with (output / "allelic_block_reassignments.tsv").open() as handle:
                events = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(events), 1)
            self.assertEqual(events[0]["acceptance_basis"], "constraint_forced")

    def test_raw_hic_mode_does_not_apply_dosage_correction(self):
        with tempfile.TemporaryDirectory() as temporary:
            output, _ = self.run_recluster(
                Path(temporary), extra=["--hic-link-normalization", "raw"]
            )
            with (output / "recluster_assignments.tsv").open() as handle:
                rows = {
                    row["unitig"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(
                rows["Y"]["decision_adjusted_assigned_hic_links"], "160.000000"
            )
            summary = json.loads((output / "recluster_summary.json").read_text())
            self.assertEqual(summary["parameters"]["hic_link_normalization"], "raw")

    def test_seed_mismatch_does_not_replace_existing_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            fasta, types, seeds, links = self.make_inputs(work)
            seeds.write_text("group1\t1\tMISSING\ngroup2\t1\tB\ngroup3\t1\tC\ngroup4\t1\tD\n")
            output = work / "output"
            output.mkdir()
            sentinel = output / "recluster_summary.json"
            sentinel.write_text("existing\n")
            result = subprocess.run(
                [
                    sys.executable,
                    str(RECLUSTER),
                    "--fasta",
                    str(fasta),
                    "--contig-type",
                    str(types),
                    "--full-links",
                    str(links),
                    "--clusters-file",
                    str(seeds),
                    "--output-dir",
                    str(output),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("absent from chromosome FASTA", result.stderr)
            self.assertEqual(sentinel.read_text(), "existing\n")

    def test_low_confidence_cluster_seed_is_reassigned_by_hic(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            fasta = work / "chr1.putg.fa"
            fasta.write_text(
                "".join(f">{name}\n" + "GATC" * 25 + "\n" for name in "ABCDQR")
            )
            types = work / "types.tsv"
            types.write_text(
                "contig_ID\taverage_depth\tcontig_type\n"
                "A\t20\thaplotig\nB\t20\thaplotig\n"
                "C\t20\thaplotig\nD\t20\thaplotig\n"
                "Q\t40\tdiplotig\nR\t20\thaplotig\n"
            )
            seeds = work / "groups.txt"
            seeds.write_text(
                "group1\t2\tA R\n"
                "group2\t2\tB Q\n"
                "group3\t2\tC Q\n"
                "group4\t1\tD\n"
            )
            cluster_assignments = work / "cluster_assignments.tsv"
            cluster_assignments.write_text(
                "unitig\tgroups\tassignment_basis\tadjusted_assigned_hic_links\t"
                "assigned_hic_fraction\thic_density_margin\n"
                "A\t1\thic_supported\t100\t1\t1\n"
                "B\t2\thic_supported\t100\t1\t1\n"
                "C\t3\thic_supported\t100\t1\t1\n"
                "D\t4\thic_supported\t100\t1\t1\n"
                "Q\t2,3\tconstraint_supported_low_margin\t2\t0.01\t-1\n"
                "R\t1\tconstraint_only_no_hic\t0\t0\t0\n"
            )
            links = work / "links.pkl"
            with links.open("wb") as handle:
                pickle.dump(
                    {
                        ("A", "Q"): 2,
                        ("B", "Q"): 2,
                        ("C", "Q"): 200,
                        ("D", "Q"): 160,
                    },
                    handle,
                )
            output = work / "output"
            subprocess.run(
                [
                    sys.executable,
                    str(RECLUSTER),
                    "--fasta", str(fasta),
                    "--contig-type", str(types),
                    "--full-links", str(links),
                    "--clusters-file", str(seeds),
                    "--cluster-assignments", str(cluster_assignments),
                    "--output-dir", str(output),
                    "--min-group-margin", "0.2",
                    "--refinement-rounds", "4",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            with (output / "recluster_assignments.tsv").open() as handle:
                rows = {row["unitig"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["Q"]["original_seed_groups"], "2,3")
            self.assertEqual(rows["Q"]["groups"], "3,4")
            self.assertEqual(rows["Q"]["seed_review_status"], "reviewed_changed")
            self.assertEqual(rows["Q"]["seed_changed"], "yes")
            self.assertEqual(
                rows["Q"]["assignment_basis"],
                "reviewed_seed_hic_high_confidence",
            )
            self.assertEqual(rows["A"]["seed_review_status"], "fixed_trusted")
            self.assertEqual(rows["R"]["groups"], "1")
            self.assertEqual(rows["R"]["seed_review_status"], "reviewed_retained")
            self.assertEqual(
                rows["R"]["assignment_basis"],
                "reviewed_seed_retained_no_assigned_hic",
            )

            summary = json.loads((output / "recluster_summary.json").read_text())
            self.assertEqual(summary["fixed_seed"]["unitigs"], 4)
            self.assertEqual(summary["reviewed_seed"]["unitigs"], 2)
            self.assertEqual(summary["reviewed_seed"]["changed_unitigs"], 1)
            self.assertEqual(summary["refinement_round_change_counts"], [0])
            self.assertTrue(summary["refinement"]["stable"])
            self.assertEqual(summary["validation"]["violations"], 0)

    def test_clm_is_split_by_common_group_in_one_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            chromosome = work / "recluster" / "chr1"
            chromosome.mkdir(parents=True)
            (chromosome / "group.reassignment.cluster.txt").write_text(
                "chr1_group1\tA X Y\n"
                "chr1_group2\tB\n"
                "chr1_group3\tC Y\n"
                "chr1_group4\tD\n"
            )
            clm = work / "links.clm"
            clm.write_text(
                "A+ X-\t10\t1 2\n"
                "Y+ C+\t8\t1 2\n"
                "Y- A+\t7\t1 2\n"
                "A+ B+\t5\t1 2\n"
            )
            subprocess.run(
                [
                    sys.executable,
                    str(SPLITTER),
                    "--clm",
                    str(clm),
                    "--recluster-dir",
                    str(work / "recluster"),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                (chromosome / "split_clms" / "group1.clm").read_text().splitlines(),
                ["A+ X-\t10\t1 2", "Y- A+\t7\t1 2"],
            )
            self.assertEqual(
                (chromosome / "split_clms" / "group3.clm").read_text().splitlines(),
                ["Y+ C+\t8\t1 2"],
            )
            summary = json.loads((work / "recluster" / "clm_split.summary.json").read_text())
            self.assertEqual(summary["input_records"], 4)
            self.assertEqual(summary["written_group_records"], 3)


if __name__ == "__main__":
    unittest.main()
