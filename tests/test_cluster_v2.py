import csv
import itertools
import json
import pickle
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLUSTER_SCRIPT = ROOT / "utils" / "cluster_allelic_unitigs_v2.py"
PUBLIC_SCRIPT = ROOT / "PHap.py"
sys.path.insert(0, str(ROOT / "utils"))
import cluster_allelic_unitigs_v2 as cluster_v2


class ConstraintClusterTests(unittest.TestCase):
    def test_public_cluster_help_has_no_clm_input(self):
        result = subprocess.run(
            [sys.executable, str(PUBLIC_SCRIPT), "cluster", "--help"],
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertNotIn("--clm", result.stdout)
        self.assertIn("--full_links", result.stdout)

    def test_hic_link_normalization_can_use_dosage_or_raw_counts(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            link_file = Path(temporary_directory) / "links.pkl"
            with link_file.open("wb") as handle:
                pickle.dump({("A", "B"): 10}, handle)

            dosage_links, _, _ = cluster_v2.load_relevant_links(
                link_file, {"A", "B"}, {"A": 1, "B": 2}, "dosage"
            )
            raw_links, _, _ = cluster_v2.load_relevant_links(
                link_file, {"A", "B"}, {"A": 1, "B": 2}, "raw"
            )

            self.assertEqual(dosage_links[("A", "B")], 5.0)
            self.assertEqual(raw_links[("A", "B")], 10.0)

    def run_cluster(
        self,
        work,
        names,
        types,
        table,
        links,
        check=True,
        extra_args=None,
        lengths=None,
    ):
        fasta = work / "chr1.fa"
        contig_types = work / "types.tsv"
        allelic_table = work / "table.tsv"
        full_links = work / "links.pkl"
        output = work / "cluster"
        lengths = lengths or {name: 100 for name in names}
        fasta.write_text(
            "".join(
                f">{name}\n"
                + (("GATC" * ((lengths[name] + 3) // 4))[: lengths[name]])
                + "\n"
                for name in names
            )
        )
        contig_types.write_text(
            "contig_ID\taverage_depth\tcontig_type\n"
            + "".join(f"{name}\t20\t{types[name]}\n" for name in names)
        )
        allelic_table.write_text(table)
        with full_links.open("wb") as handle:
            pickle.dump(links, handle)
        command = [
                sys.executable,
                str(CLUSTER_SCRIPT),
                "--fasta",
                str(fasta),
                "--full-links",
                str(full_links),
                "--contig-type",
                str(contig_types),
                "--allelic-table",
                str(allelic_table),
                "--output-dir",
                str(output),
            ]
        if extra_args:
            command.extend(extra_args)
        result = subprocess.run(
            command,
            check=check,
            capture_output=True,
            text=True,
        )
        return output, result

    def test_enforces_allelic_constraints_and_dosage_without_dropping_no_hic_unitigs(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            names = ["A", "B", "C", "D", "E", "X", "Y", "Z", "N"]
            types = {name: "haplotig" for name in names}
            types["X"] = "diplotig"
            table = (
                "chr1\t0\t100\tA\tB\tC\tD\n"
                "chr1\t100\t200\tA\tE\n"
                "chr1\t200\t300\tX\tY\tZ\n"
                "chr1\t300\t400\tN\n"
            )
            links = {
                tuple(sorted(("E", "B"))): 1000,
                tuple(sorted(("E", "C"))): 10,
                tuple(sorted(("X", "A"))): 500,
                tuple(sorted(("X", "B"))): 500,
                tuple(sorted(("Y", "C"))): 500,
                tuple(sorted(("Z", "D"))): 500,
            }
            output, result = self.run_cluster(work, names, types, table, links)
            self.assertIn("clustered 9 table unitigs", result.stdout)

            with (output / "cluster_assignments.tsv").open() as handle:
                rows = {row["unitig"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(set(rows), set(names))
            memberships = {
                unitig: set(map(int, row["groups"].split(",")))
                for unitig, row in rows.items()
            }
            self.assertEqual(len(set(next(iter(memberships[item])) for item in "ABCD")), 4)
            self.assertTrue(memberships["A"].isdisjoint(memberships["E"]))
            self.assertEqual(memberships["E"], memberships["B"])
            self.assertEqual(len(memberships["X"]), 2)
            self.assertTrue(memberships["X"].isdisjoint(memberships["Y"]))
            self.assertTrue(memberships["X"].isdisjoint(memberships["Z"]))
            self.assertTrue(memberships["Y"].isdisjoint(memberships["Z"]))
            self.assertEqual(rows["N"]["assignment_basis"], "constraint_only_no_hic")

            summary = json.loads((output / "cluster_summary.json").read_text())
            self.assertEqual(summary["parameters"]["hic_link_normalization"], "dosage")
            self.assertEqual(summary["validation"]["allelic_conflicts"], 0)
            self.assertEqual(summary["validation"]["dosage_errors"], 0)
            self.assertEqual(summary["validation"]["unassigned_table_unitigs"], 0)
            self.assertGreaterEqual(
                summary["refinement"]["objective_final"],
                summary["refinement"]["objective_initial"],
            )
            self.assertIn("kempe_moves", summary["refinement"])
            self.assertIn("phase_block_moves", summary["refinement"])
            self.assertTrue((output / "cluster_phase_block_moves.tsv").exists())
            self.assertTrue(
                (output / "cluster_constraint_relaxation_moves.tsv").exists()
            )
            self.assertTrue((output / "cluster_phase_interval_moves.tsv").exists())
            self.assertEqual(
                (output / "cluster_constraint_violations.tsv").read_text().splitlines(),
                ["violation\tunitig1\tunitig2\texpected\tobserved"],
            )
            self.assertFalse((output / "full.links.txt").exists())
            self.assertFalse((output / "flank.links.txt").exists())

    def test_relaxes_weakest_edge_in_globally_overcapacity_clique(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            names = ["A", "B", "C", "D", "E"]
            types = {name: "haplotig" for name in names}
            table_rows = []
            coordinate = 0
            for index, unitig1 in enumerate(names):
                for unitig2 in names[index + 1:]:
                    table_rows.append(
                        f"chr1\t{coordinate}\t{coordinate + 10}\t{unitig1}\t{unitig2}\n"
                    )
                    coordinate += 10
            output, result = self.run_cluster(
                work,
                names,
                types,
                "".join(table_rows),
                {},
            )
            self.assertEqual(result.returncode, 0)
            summary = json.loads((output / "cluster_summary.json").read_text())
            self.assertEqual(summary["table"]["relaxed_conflict_pairs"], 1)
            self.assertEqual(summary["validation"]["allelic_conflicts"], 0)
            self.assertEqual(summary["validation"]["dosage_errors"], 0)
            relaxed = (output / "cluster_relaxed_constraints.tsv").read_text().splitlines()
            self.assertEqual(len(relaxed), 2)

    def test_relaxation_strength_uses_sequence_length_not_repeated_table_bp(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            names = ["A", "B", "C", "D", "E"]
            types = {name: "haplotig" for name in names}
            table_rows = ["chr1\t0\t5\tA\tB\n"]
            coordinate = 10
            for left, right in itertools.combinations(names, 2):
                if (left, right) == ("A", "B"):
                    continue
                table_rows.append(
                    f"chr1\t{coordinate}\t{coordinate + 10}\t{left}\t{right}\n"
                )
                coordinate += 10
            table_rows.append("chr1\t1000\t10000\tE\n")

            output, result = self.run_cluster(
                work,
                names,
                types,
                "".join(table_rows),
                {},
                lengths={name: 100 for name in names},
            )

            self.assertEqual(result.returncode, 0)
            with (output / "cluster_relaxed_constraints.tsv").open() as handle:
                relaxed = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                {relaxed[0]["unitig1"], relaxed[0]["unitig2"]}, {"A", "B"}
            )
            self.assertAlmostEqual(float(relaxed[0]["ratio1"]), 0.05)
            self.assertAlmostEqual(float(relaxed[0]["ratio2"]), 0.05)

    def test_protects_long_unitig_containment_edges_during_relaxation(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            names = ["LONG", "B", "C", "D", "E"]
            types = {name: "haplotig" for name in names}
            table = "".join(
                f"chr1\t{index * 100}\t{index * 100 + 100}\t{left}\t{right}\n"
                for index, (left, right) in enumerate(
                    itertools.combinations(names, 2)
                )
            )
            output, result = self.run_cluster(
                work,
                names,
                types,
                table,
                {},
                lengths={"LONG": 1000, "B": 100, "C": 100, "D": 100, "E": 100},
                extra_args=[
                    "--protected-long-unitig-length",
                    "500",
                    "--protected-length-ratio",
                    "5",
                    "--protected-short-overlap",
                    "0.5",
                ],
            )

            self.assertEqual(result.returncode, 0)
            summary = json.loads((output / "cluster_summary.json").read_text())
            self.assertEqual(summary["table"]["protected_conflict_pairs"], 4)
            with (output / "cluster_relaxed_constraints.tsv").open() as handle:
                relaxed = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(relaxed)
            self.assertTrue(
                all("LONG" not in (row["unitig1"], row["unitig2"]) for row in relaxed)
            )
            with (output / "cluster_protected_constraints.tsv").open() as handle:
                protected = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(protected), 4)
            self.assertTrue(all(row["long_unitig"] == "LONG" for row in protected))

    def test_strict_mode_rejects_globally_infeasible_constraints(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            names = ["A", "B", "C", "D", "E"]
            types = {name: "haplotig" for name in names}
            table = "".join(
                f"chr1\t{index * 10}\t{index * 10 + 10}\t{unitig1}\t{unitig2}\n"
                for index, (unitig1, unitig2) in enumerate(
                    itertools.combinations(names, 2)
                )
            )
            output, result = self.run_cluster(
                work,
                names,
                types,
                table,
                {},
                check=False,
                extra_args=["--constraint-relaxation", "fail"],
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("not satisfiable", result.stderr)
            self.assertFalse((output / "cluster_summary.json").exists())

    def test_relaxes_non_clique_constraint_component_that_is_not_colorable(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            work = Path(temporary_directory)
            names = ["A", "B", "C", "D", "E"]
            types = {name: "haplotig" for name in names}
            cycle = [("A", "B"), ("B", "C"), ("C", "D"), ("D", "E"), ("E", "A")]
            table = "".join(
                f"chr1\t{index * 10}\t{index * 10 + 10}\t{left}\t{right}\n"
                for index, (left, right) in enumerate(cycle)
            )
            output, result = self.run_cluster(
                work,
                names,
                types,
                table,
                {},
                extra_args=["--ploidy", "2"],
            )

            self.assertEqual(result.returncode, 0)
            summary = json.loads((output / "cluster_summary.json").read_text())
            self.assertEqual(summary["table"]["relaxed_conflict_pairs"], 1)
            with (output / "cluster_relaxed_constraints.tsv").open() as handle:
                relaxed = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                relaxed[0]["reason"], "constraint_component_unsatisfiable"
            )

    def test_phase_block_refinement_swaps_suffix_and_relaxes_one_weak_boundary(self):
        units = ["A", "B", "C", "D", "E", "F"]
        dosage = {unitig: 1 for unitig in units}
        lengths = {unitig: 100 for unitig in units}
        re_sites = {unitig: 10 for unitig in units}
        adjacency = {unitig: set() for unitig in units}
        adjacency["A"].add("E")
        adjacency["E"].update(("A", "F"))
        adjacency["F"].add("E")
        links = {("A", "E"): 1000.0, ("D", "F"): 1000.0}
        link_neighbors = {
            "A": {"E": 1000.0},
            "E": {"A": 1000.0},
            "D": {"F": 1000.0},
            "F": {"D": 1000.0},
        }
        rows = [
            cluster_v2.TableRow("chr1", 0, 10, ("A", "B", "C", "D")),
            cluster_v2.TableRow("chr1", 40, 41, ("A", "E")),
            cluster_v2.TableRow("chr1", 90, 100, ("E", "F")),
        ]
        clusterer = cluster_v2.ConstraintClusterer(
            units,
            dosage,
            lengths,
            re_sites,
            adjacency,
            links,
            link_neighbors,
            rows,
            4,
            1.0,
            1000,
        )
        assignments = {
            "A": 1 << 0,
            "B": 1 << 1,
            "C": 1 << 2,
            "D": 1 << 3,
            "E": 1 << 3,
            "F": 1 << 0,
        }
        edge_overlap_bp = {("A", "E"): 1, ("E", "F"): 10}

        assignments, metrics, relaxed = clusterer.refine_phase_blocks(
            assignments,
            max_rounds=1,
            min_objective_gain=0.005,
            max_boundary_relaxations=1,
            max_boundary_overlap=0.30,
            edge_overlap_bp=edge_overlap_bp,
            protected_edges=set(),
        )

        self.assertEqual(assignments["E"], 1 << 0)
        self.assertEqual(assignments["F"], 1 << 3)
        self.assertEqual(len(metrics["moves"]), 1)
        self.assertGreater(metrics["final_score"], metrics["initial_score"])
        self.assertEqual(
            [(item["unitig1"], item["unitig2"], item["reason"]) for item in relaxed],
            [("A", "E", "phase_block_boundary")],
        )
        self.assertNotIn("E", adjacency["A"])
        self.assertEqual(
            cluster_v2.validate_assignments(assignments, dosage, adjacency, 4), []
        )

    def test_strong_hic_can_release_only_a_weak_unprotected_constraint(self):
        units = ["A", "B", "X"]
        dosage = {unitig: 1 for unitig in units}
        lengths = {unitig: 100 for unitig in units}
        re_sites = {unitig: 10 for unitig in units}
        adjacency = {unitig: set() for unitig in units}
        adjacency["A"].add("X")
        adjacency["X"].add("A")
        links = {("A", "X"): 1000.0}
        link_neighbors = {
            "A": {"X": 1000.0},
            "X": {"A": 1000.0},
        }
        clusterer = cluster_v2.ConstraintClusterer(
            units,
            dosage,
            lengths,
            re_sites,
            adjacency,
            links,
            link_neighbors,
            [cluster_v2.TableRow("chr1", 0, 10, ("A", "X"))],
            2,
            1.0,
            1000,
        )
        assignments = {"A": 1 << 0, "B": 1 << 1, "X": 1 << 1}

        assignments, metrics, relaxed = clusterer.refine_contradicted_constraints(
            assignments,
            max_rounds=2,
            min_objective_gain=0.0,
            min_adjusted_links=5.0,
            min_group_margin=0.10,
            max_relaxations=1,
            max_overlap=0.30,
            edge_overlap_bp={("A", "X"): 10},
            protected_edges=set(),
        )

        self.assertEqual(assignments["X"], 1 << 0)
        self.assertEqual(len(metrics["moves"]), 1)
        self.assertEqual(len(relaxed), 1)
        self.assertEqual(relaxed[0]["reason"], "hic_contradicted_weak_constraint")
        self.assertNotIn("X", adjacency["A"])

    def test_phase_interval_refinement_swaps_an_internal_block(self):
        units = ["A", "B", "E", "F", "C", "D"]
        dosage = {unitig: 1 for unitig in units}
        lengths = {unitig: 100 for unitig in units}
        re_sites = {unitig: 10 for unitig in units}
        adjacency = {unitig: set() for unitig in units}
        for left, right in (("A", "B"), ("E", "F"), ("C", "D")):
            adjacency[left].add(right)
            adjacency[right].add(left)
        links = {
            ("A", "E"): 100.0,
            ("B", "F"): 100.0,
            ("C", "E"): 100.0,
            ("D", "F"): 100.0,
        }
        link_neighbors = {
            "A": {"E": 100.0},
            "B": {"F": 100.0},
            "C": {"E": 100.0},
            "D": {"F": 100.0},
            "E": {"A": 100.0, "C": 100.0},
            "F": {"B": 100.0, "D": 100.0},
        }
        rows = [
            cluster_v2.TableRow("chr1", 0, 10, ("A", "B")),
            cluster_v2.TableRow("chr1", 50, 60, ("E", "F")),
            cluster_v2.TableRow("chr1", 100, 110, ("C", "D")),
        ]
        clusterer = cluster_v2.ConstraintClusterer(
            units,
            dosage,
            lengths,
            re_sites,
            adjacency,
            links,
            link_neighbors,
            rows,
            2,
            1.0,
            1000,
        )
        assignments = {
            "A": 1 << 0,
            "B": 1 << 1,
            "E": 1 << 1,
            "F": 1 << 0,
            "C": 1 << 0,
            "D": 1 << 1,
        }
        edge_overlap_bp = {
            ("A", "B"): 10,
            ("E", "F"): 10,
            ("C", "D"): 10,
        }

        assignments, metrics, relaxed = clusterer.refine_phase_intervals(
            assignments,
            max_rounds=1,
            min_objective_gain=0.005,
            max_boundary_relaxations=0,
            max_boundary_overlap=0.30,
            edge_overlap_bp=edge_overlap_bp,
            protected_edges=set(),
        )

        self.assertEqual(assignments["E"], 1 << 0)
        self.assertEqual(assignments["F"], 1 << 1)
        self.assertEqual(len(metrics["moves"]), 1)
        self.assertEqual(relaxed, [])
        self.assertGreater(metrics["final_score"], metrics["initial_score"])
        self.assertEqual(
            cluster_v2.validate_assignments(assignments, dosage, adjacency, 2), []
        )


if __name__ == "__main__":
    unittest.main()
