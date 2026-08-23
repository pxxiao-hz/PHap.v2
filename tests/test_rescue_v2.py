import csv
import json
import pickle
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESCUE = ROOT / "utils" / "unchr_recluster.py"
SPLITTER = ROOT / "utils" / "split_clm_by_groups_v2.py"


class RescueTests(unittest.TestCase):
    def write_chromosome(self, root, chromosome, seeds, deferred=None):
        deferred = deferred or []
        directory = root / chromosome
        directory.mkdir(parents=True)
        header = (
            "unitig\tlength\tre_sites\tcontig_type\tdosage\tdosage_source\t"
            "status\tgroups\tassignment_basis\n"
        )
        rows = []
        for group, unitig in enumerate(seeds, 1):
            rows.append(
                f"{unitig}\t20\t6\thaplotig\t1\tprovided\tassigned\t"
                f"{group}\tfixed_seed\n"
            )
            (directory / f"group{group}.reassignment.fa").write_text(
                f">{unitig}\n" + "GATC" * 5 + "\n"
            )
        for unitig in deferred:
            rows.append(
                f"{unitig}\t20\t6\thaplotig\t1\tprovided\tunassigned\t"
                "\tno_assigned_hic\n"
            )
        (directory / "recluster_assignments.tsv").write_text(header + "".join(rows))
        (directory / "group.reassignment.cluster.txt").write_text(
            "".join(
                f"{chromosome}_group{group}\t{unitig}\n"
                for group, unitig in enumerate(seeds, 1)
            )
        )
        (directory / "unassigned_unitigs.fa").write_text(
            "".join(f">{unitig}\n" + "GATC" * 5 + "\n" for unitig in deferred)
        )

    def make_inputs(self, work):
        recluster = work / "04.recluster"
        self.write_chromosome(recluster, "chr01", ["A", "B", "C", "D"], ["Z"])
        self.write_chromosome(recluster, "chr02", ["E", "F", "G", "H"])
        candidates = ["X", "Y", "R", "P", "O"]
        candidate_fasta = work / "un_chr.fa"
        candidate_fasta.write_text(
            "".join(f">{unitig}\n" + "GATC" * 5 + "\n" for unitig in candidates)
        )
        assembly = work / "assembly.fa"
        all_unitigs = ["A", "B", "C", "D", "Z", "E", "F", "G", "H"] + candidates
        assembly.write_text(
            "".join(f">{unitig}\n" + "GATC" * 5 + "\n" for unitig in all_unitigs)
        )
        types = work / "types.tsv"
        candidate_types = {
            "X": "haplotig",
            "Y": "diplotig",
            "R": "haplotig",
            "P": "haplotig",
            "O": "other",
        }
        types.write_text(
            "contig_ID\taverage_depth\tcontig_type\n"
            + "".join(
                f"{unitig}\t20\t{candidate_types[unitig]}\n" for unitig in candidates
            )
        )
        links = {
            tuple(sorted(("X", "A"))): 100,
            tuple(sorted(("R", "X"))): 100,
            tuple(sorted(("Y", "F"))): 100,
            tuple(sorted(("Y", "G"))): 100,
            tuple(sorted(("P", "A"))): 50,
            tuple(sorted(("P", "E"))): 50,
            tuple(sorted(("P", "X"))): 50,
            tuple(sorted(("P", "R"))): 50,
            tuple(sorted(("O", "A"))): 100,
            tuple(sorted(("Z", "A"))): 100,
        }
        link_file = work / "links.pkl"
        with link_file.open("wb") as handle:
            pickle.dump(links, handle)
        return assembly, candidate_fasta, types, link_file, recluster

    def run_rescue(self, work, extra=None, check=True):
        assembly, candidates, types, links, recluster = self.make_inputs(work)
        output = work / "05.rescue"
        command = [
            sys.executable,
            str(RESCUE),
            "--assembly-fasta", str(assembly),
            "--candidate-fasta", str(candidates),
            "--contig-type", str(types),
            "--full-links", str(links),
            "--recluster-dir", str(recluster),
            "--output-dir", str(output),
            "--min-chromosome-margin", "0.2",
            "--min-group-margin", "0.2",
        ]
        if extra:
            command.extend(extra)
        result = subprocess.run(command, check=check, capture_output=True, text=True)
        return output, result

    def test_candidate_scope_exact_dosage_and_propagation(self):
        with tempfile.TemporaryDirectory() as temporary:
            output, result = self.run_rescue(Path(temporary))
            self.assertIn("Rescued 3/5", result.stdout)
            with (output / "rescue_assignments.tsv").open() as handle:
                rows = {
                    row["unitig"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }

            self.assertEqual(rows["X"]["groups"], "chr01_group1")
            self.assertEqual(rows["R"]["groups"], "chr01_group1")
            self.assertEqual(rows["R"]["assignment_round"], "2")
            self.assertEqual(rows["Y"]["groups"], "chr02_group2,chr02_group3")
            self.assertEqual(rows["Y"]["group_count"], rows["Y"]["dosage"])
            self.assertEqual(
                rows["Y"]["decision_adjusted_selected_hic_links"], "100.000000"
            )
            self.assertEqual(rows["P"]["assignment_basis"], "low_chromosome_margin")
            self.assertEqual(rows["O"]["assignment_basis"], "unsupported_contig_type:other")
            self.assertEqual(rows["Z"]["status"], "unassigned")
            self.assertEqual(rows["Z"]["source"], "recluster_deferred")

            summary = json.loads((output / "rescue_summary.json").read_text())
            self.assertEqual(summary["parameters"]["hic_link_normalization"], "dosage")
            self.assertEqual(summary["candidates"]["rescued_unitigs"], 3)
            self.assertEqual(summary["round_assignment_counts"], [2, 1])
            self.assertEqual(summary["validation"]["violations"], 0)
            self.assertEqual(summary["final"]["assigned_unitigs"], 11)
            self.assertEqual(
                set((output / "unassigned_unitigs.txt").read_text().split()),
                {"Z", "P", "O"},
            )
            self.assertIn(">X\n", (output / "chr01_group1.reassignment.fa").read_text())
            self.assertFalse((output / "g1g2g3g4.reassignment.fa").exists())

    def test_raw_hic_mode_does_not_apply_dosage_correction(self):
        with tempfile.TemporaryDirectory() as temporary:
            output, _ = self.run_rescue(
                Path(temporary), extra=["--hic-link-normalization", "raw"]
            )
            with (output / "rescue_assignments.tsv").open() as handle:
                rows = {
                    row["unitig"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(
                rows["Y"]["decision_adjusted_selected_hic_links"], "200.000000"
            )
            summary = json.loads((output / "rescue_summary.json").read_text())
            self.assertEqual(summary["parameters"]["hic_link_normalization"], "raw")

    def test_input_failure_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            assembly, candidates, types, links, recluster = self.make_inputs(work)
            (recluster / "chr01" / "group1.reassignment.fa").write_text(
                ">WRONG\nGATC\n"
            )
            output = work / "05.rescue"
            output.mkdir()
            sentinel = output / "rescue_summary.json"
            sentinel.write_text("existing\n")
            result = subprocess.run(
                [
                    sys.executable, str(RESCUE),
                    "--assembly-fasta", str(assembly),
                    "--candidate-fasta", str(candidates),
                    "--contig-type", str(types),
                    "--full-links", str(links),
                    "--recluster-dir", str(recluster),
                    "--output-dir", str(output),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("FASTA membership/order mismatch", result.stderr)
            self.assertEqual(sentinel.read_text(), "existing\n")

    def test_final_cluster_file_clm_split(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            output, _ = self.run_rescue(work)
            clm = work / "links.clm"
            clm.write_text(
                "A+ X-\t10\t1 2\n"
                "F+ Y+\t8\t1 2\n"
                "G- Y+\t7\t1 2\n"
                "A+ E+\t5\t1 2\n"
            )
            subprocess.run(
                [
                    sys.executable, str(SPLITTER),
                    "--clm", str(clm),
                    "--clusters-file", str(output / "group.reassignment.cluster.txt"),
                    "--output-dir", str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                (output / "split_clms" / "chr01_group1.clm").read_text().splitlines(),
                ["A+ X-\t10\t1 2"],
            )
            self.assertEqual(
                (output / "split_clms" / "chr02_group2.clm").read_text().splitlines(),
                ["F+ Y+\t8\t1 2"],
            )
            summary = json.loads((output / "clm_split.summary.json").read_text())
            self.assertEqual(summary["input_records"], 4)
            self.assertEqual(summary["written_group_records"], 3)


if __name__ == "__main__":
    unittest.main()
