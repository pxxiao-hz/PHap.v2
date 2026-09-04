import csv
import importlib.util
import json
import pickle
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
UTILS = ROOT / "utils"
sys.path.insert(0, str(UTILS))

from dosage import contig_type_for_dosage, dosage_from_contig_type


class GeneralPloidyTests(unittest.TestCase):
    def test_dosage_labels_are_backward_compatible_and_extensible(self):
        expected = {
            "haplotig": 1,
            "triplotig": 3,
            "pentaplotig": 5,
            "hexaplotig": 6,
            "dosage_7": 7,
            "copy-8": 8,
            "9x": 9,
            "10": 10,
        }
        for label, dosage in expected.items():
            with self.subTest(label=label):
                self.assertEqual(dosage_from_contig_type(label), dosage)
        self.assertIsNone(dosage_from_contig_type("replotig"))
        self.assertEqual(contig_type_for_dosage(5), "pentaplotig")
        self.assertEqual(contig_type_for_dosage(7), "dosage_7")

    def test_depth_classifier_uses_ploidy_and_single_copy_depth(self):
        path = ROOT / "shell" / "dosage.analysis.contig.type.identified.py"
        spec = importlib.util.spec_from_file_location("depth_classifier", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        self.assertEqual(module.classify_contig_normal(30, 20, 3), "diplotig")
        self.assertEqual(module.classify_contig_normal(60, 20, 3), "triplotig")
        self.assertEqual(module.classify_contig_normal(100, 20, 6), "pentaplotig")
        self.assertEqual(module.classify_contig_normal(120, 20, 6), "hexaplotig")
        self.assertEqual(module.classify_contig_normal(140, 20, 6), "replotig")

    def test_no_collapse_cluster_needs_no_contig_type_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            names = ["A", "B", "C", "X"]
            fasta = work / "chr1.fa"
            fasta.write_text(
                "".join(f">{name}\n" + "GATC" * 5 + "\n" for name in names)
            )
            table = work / "table.tsv"
            table.write_text(
                "chr1\t0\t100\tA\tB\tC\n"
                "chr1\t100\t200\tX\n"
            )
            links = work / "links.pkl"
            with links.open("wb") as handle:
                pickle.dump({("A", "X"): 100}, handle)
            output = work / "cluster"
            subprocess.run(
                [
                    sys.executable,
                    str(UTILS / "cluster_allelic_unitigs_v2.py"),
                    "--fasta", str(fasta),
                    "--full-links", str(links),
                    "--allelic-table", str(table),
                    "--output-dir", str(output),
                    "--ploidy", "3",
                    "--no-collapse",
                ],
                check=True,
            )
            with (output / "cluster_assignments.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(all(row["dosage"] == "1" for row in rows))
            self.assertTrue(all(row["group_count"] == "1" for row in rows))
            summary = json.loads((output / "cluster_summary.json").read_text())
            self.assertTrue(summary["parameters"]["no_collapse"])
            self.assertIsNone(summary["inputs"]["contig_type"])

    def test_no_collapse_allelic_table_needs_no_contig_type_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            paf = work / "selected.paf"
            paf.write_text(
                "A\t100000\t0\t100000\t+\tchr1\t1000000\t0\t100000\t"
                "98000\t100000\t60\ttp:A:P\n"
            )
            paths = {
                "output": work / "table.tsv",
                "projections": work / "projections.tsv",
                "qc": work / "qc.tsv",
                "rejected": work / "rejected.tsv",
                "pairs": work / "pairs.tsv",
                "summary": work / "summary.json",
            }
            command = [
                sys.executable,
                str(UTILS / "allelic_table_generate_v2.py"),
                "--paf", str(paf),
                "--ploidy", "3",
                "--no-collapse",
            ]
            for option, path in paths.items():
                command.extend([f"--{option}", str(path)])
            subprocess.run(command, check=True)
            with paths["projections"].open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["contig_type"], "haplotig")
            self.assertEqual(row["dosage"], "1")
            summary = json.loads(paths["summary"].read_text())
            self.assertTrue(summary["parameters"]["no_collapse"])

    def test_hexaploid_cluster_assigns_exact_high_dosages(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            names = list("ABCDEF") + ["P", "H"]
            fasta = work / "chr1.fa"
            fasta.write_text("".join(f">{name}\n" + "GATC" * 25 + "\n" for name in names))
            types = work / "types.tsv"
            types.write_text(
                "contig_ID\taverage_depth\tcontig_type\n"
                + "".join(f"{name}\t20\thaplotig\n" for name in "ABCDEF")
                + "P\t100\tpentaplotig\n"
                + "H\t120\tdosage_6\n"
            )
            table = work / "table.tsv"
            table.write_text(
                "chr1\t0\t100\tA\tB\tC\tD\tE\tF\n"
                "chr1\t100\t200\tP\n"
                "chr1\t200\t300\tH\n"
            )
            links = work / "links.pkl"
            with links.open("wb") as handle:
                pickle.dump({}, handle)
            output = work / "cluster"
            subprocess.run(
                [
                    sys.executable,
                    str(UTILS / "cluster_allelic_unitigs_v2.py"),
                    "--fasta", str(fasta),
                    "--full-links", str(links),
                    "--contig-type", str(types),
                    "--allelic-table", str(table),
                    "--output-dir", str(output),
                    "--ploidy", "6",
                ],
                check=True,
            )
            with (output / "cluster_assignments.tsv").open() as handle:
                rows = {row["unitig"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["P"]["dosage"], "5")
            self.assertEqual(rows["P"]["group_count"], "5")
            self.assertEqual(rows["H"]["dosage"], "6")
            self.assertEqual(rows["H"]["group_count"], "6")
            self.assertTrue((output / "g6.fa").exists())
            self.assertTrue((output / "all_groups.fa").exists())
            self.assertFalse((output / "g1g2g3g4.fa").exists())
            summary = json.loads((output / "cluster_summary.json").read_text())
            self.assertEqual(summary["parameters"]["ploidy"], 6)
            self.assertEqual(summary["validation"]["dosage_errors"], 0)

    def test_triploid_and_pentaploid_cluster_group_counts(self):
        for ploidy, full_type in ((3, "triplotig"), (5, "pentaplotig")):
            with self.subTest(ploidy=ploidy), tempfile.TemporaryDirectory() as temporary:
                work = Path(temporary)
                seeds = [f"S{number}" for number in range(1, ploidy + 1)]
                names = seeds + ["FULL"]
                fasta = work / "chr1.fa"
                fasta.write_text(
                    "".join(f">{name}\n" + "GATC" * 5 + "\n" for name in names)
                )
                types = work / "types.tsv"
                types.write_text(
                    "contig_ID\taverage_depth\tcontig_type\n"
                    + "".join(f"{name}\t20\thaplotig\n" for name in seeds)
                    + f"FULL\t{20 * ploidy}\t{full_type}\n"
                )
                table = work / "table.tsv"
                seed_fields = "\t".join(seeds)
                table.write_text(
                    f"chr1\t0\t100\t{seed_fields}\n"
                    "chr1\t100\t200\tFULL\n"
                )
                links = work / "links.pkl"
                with links.open("wb") as handle:
                    pickle.dump({}, handle)
                output = work / "cluster"
                subprocess.run(
                    [
                        sys.executable,
                        str(UTILS / "cluster_allelic_unitigs_v2.py"),
                        "--fasta", str(fasta),
                        "--full-links", str(links),
                        "--contig-type", str(types),
                        "--allelic-table", str(table),
                        "--output-dir", str(output),
                        "--ploidy", str(ploidy),
                    ],
                    check=True,
                )
                with (output / "cluster_assignments.tsv").open() as handle:
                    rows = {
                        row["unitig"]: row
                        for row in csv.DictReader(handle, delimiter="\t")
                    }
                self.assertEqual(rows["FULL"]["dosage"], str(ploidy))
                self.assertEqual(rows["FULL"]["group_count"], str(ploidy))
                self.assertTrue((output / f"g{ploidy}.fa").exists())

    def test_allelic_table_accepts_pentaplotig_at_ploidy_six(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            paf = work / "selected.paf"
            paf.write_text(
                "P\t100000\t0\t100000\t+\tchr1\t1000000\t0\t100000\t"
                "98000\t100000\t60\ttp:A:P\n"
            )
            types = work / "types.tsv"
            types.write_text("contig_ID\tcontig_type\nP\tpentaplotig\n")
            paths = {
                "output": work / "table.tsv",
                "projections": work / "projections.tsv",
                "qc": work / "qc.tsv",
                "rejected": work / "rejected.tsv",
                "pairs": work / "pairs.tsv",
                "summary": work / "summary.json",
            }
            command = [
                sys.executable,
                str(UTILS / "allelic_table_generate_v2.py"),
                "--paf", str(paf),
                "--contig-type", str(types),
                "--ploidy", "6",
            ]
            for option, path in paths.items():
                command.extend([f"--{option}", str(path)])
            subprocess.run(command, check=True)
            with paths["projections"].open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["contig_type"], "pentaplotig")
            self.assertEqual(row["dosage"], "5")
            self.assertIn("\tP\n", paths["output"].read_text())

    def test_hexaploid_recluster_assigns_full_dosage_without_hic(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            names = list("ABCDEF") + ["H"]
            fasta = work / "chr1.putg.fa"
            fasta.write_text("".join(f">{name}\n" + "GATC" * 5 + "\n" for name in names))
            types = work / "types.tsv"
            types.write_text(
                "contig_ID\taverage_depth\tcontig_type\n"
                + "".join(f"{name}\t20\thaplotig\n" for name in "ABCDEF")
                + "H\t120\thexaplotig\n"
            )
            seeds = work / "groups.txt"
            seeds.write_text(
                "".join(
                    f"group{number}\t1\t{name}\n"
                    for number, name in enumerate("ABCDEF", 1)
                )
            )
            links = work / "links.pkl"
            with links.open("wb") as handle:
                pickle.dump({}, handle)
            output = work / "recluster"
            subprocess.run(
                [
                    sys.executable,
                    str(UTILS / "chr_uncluster_recluster.py"),
                    "--fasta", str(fasta),
                    "--contig-type", str(types),
                    "--full-links", str(links),
                    "--clusters-file", str(seeds),
                    "--output-dir", str(output),
                    "--ploidy", "6",
                    "--min-group-bp-ratio", "0",
                ],
                check=True,
            )
            with (output / "recluster_assignments.tsv").open() as handle:
                rows = {row["unitig"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["H"]["groups"], "1,2,3,4,5,6")
            self.assertEqual(rows["H"]["assignment_basis"], "dosage_all_groups")
            self.assertTrue((output / "group6.reassignment.fa").exists())
            self.assertTrue((output / "all_groups.reassignment.fa").exists())

    def test_hexaploid_rescue_uses_all_six_groups(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            recluster = work / "04.recluster" / "chr1"
            recluster.mkdir(parents=True)
            header = (
                "unitig\tlength\tre_sites\tcontig_type\tdosage\tdosage_source\t"
                "status\tgroups\tassignment_basis\n"
            )
            rows = []
            group_lines = []
            for number, name in enumerate("ABCDEF", 1):
                rows.append(
                    f"{name}\t20\t6\thaplotig\t1\tprovided\tassigned\t"
                    f"{number}\tfixed_seed\n"
                )
                group_lines.append(f"chr1_group{number}\t{name}\n")
                (recluster / f"group{number}.reassignment.fa").write_text(
                    f">{name}\n" + "GATC" * 5 + "\n"
                )
            (recluster / "recluster_assignments.tsv").write_text(header + "".join(rows))
            (recluster / "group.reassignment.cluster.txt").write_text("".join(group_lines))
            (recluster / "unassigned_unitigs.fa").write_text("")

            candidates = work / "un_chr.fa"
            candidates.write_text(">H\n" + "GATC" * 5 + "\n")
            assembly = work / "assembly.fa"
            assembly.write_text(
                "".join(f">{name}\n" + "GATC" * 5 + "\n" for name in "ABCDEFH")
            )
            types = work / "types.tsv"
            types.write_text(
                "contig_ID\taverage_depth\tcontig_type\nH\t120\thexaplotig\n"
            )
            links = work / "links.pkl"
            with links.open("wb") as handle:
                pickle.dump({tuple(sorted(("H", name))): 100 for name in "ABCDEF"}, handle)
            output = work / "05.rescue"
            subprocess.run(
                [
                    sys.executable,
                    str(UTILS / "unchr_recluster.py"),
                    "--assembly-fasta", str(assembly),
                    "--candidate-fasta", str(candidates),
                    "--contig-type", str(types),
                    "--full-links", str(links),
                    "--recluster-dir", str(recluster.parent),
                    "--output-dir", str(output),
                    "--ploidy", "6",
                ],
                check=True,
            )
            with (output / "rescue_assignments.tsv").open() as handle:
                result = {row["unitig"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(result["H"]["dosage"], "6")
            self.assertEqual(result["H"]["group_count"], "6")
            self.assertEqual(
                result["H"]["groups"],
                "chr1_group1,chr1_group2,chr1_group3,chr1_group4,chr1_group5,chr1_group6",
            )

    def test_phase_reads_validates_six_group_membership(self):
        from phase_reads_assignment import parse_contig_types, parse_group_file, validate_group_dosage

        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            groups = work / "groups.txt"
            groups.write_text(
                "".join(f"chr1_group{number}\tH\n" for number in range(1, 7))
            )
            types = work / "types.tsv"
            types.write_text(
                "contig_ID\taverage_depth\tcontig_type\nH\t120\thexaplotig\n"
            )
            model = parse_group_file(groups)
            summary = validate_group_dosage(model, parse_contig_types(types))
            self.assertEqual(model.ploidy, 6)
            self.assertEqual(summary["validated"], 1)


if __name__ == "__main__":
    unittest.main()
