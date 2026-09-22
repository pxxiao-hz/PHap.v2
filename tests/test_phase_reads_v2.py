import gzip
import json
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace

import pysam


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "utils"))
import phase_reads_assignment as phase
import phase_reads_assemble_anchor as phase_workflow

PHASE_SCRIPT = ROOT / "utils" / "phase_reads_assemble_anchor.py"


class PhaseReadsAssignmentTests(unittest.TestCase):
    def make_model(self):
        return phase.GroupModel(
            groups={
                "chr1_group1": ("A", "C"),
                "chr1_group2": ("B", "C"),
            },
            unitig_groups={
                "A": ("chr1_group1",),
                "B": ("chr1_group2",),
                "C": ("chr1_group1", "chr1_group2"),
            },
            group_chromosome={
                "chr1_group1": "chr1",
                "chr1_group2": "chr1",
            },
            group_number={"chr1_group1": 1, "chr1_group2": 2},
            ploidy=2,
        )

    def hit(self, unitig, groups, length=1000):
        return phase.AlignmentHit(unitig, 0, length, 1.0, 60, tuple(groups))

    def test_group_parser_and_dosage_validation(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            groups = work / "groups.txt"
            groups.write_text(
                "chr1_group1\tA C\n"
                "chr1_group2\tB C\n"
            )
            types = work / "types.tsv"
            types.write_text(
                "contig_ID\taverage_depth\tcontig_type\n"
                "A\t20\thaplotig\n"
                "B\t20\thaplotig\n"
                "C\t40\tdiplotig\n"
            )
            model = phase.parse_group_file(groups)
            summary = phase.validate_group_dosage(
                model, phase.parse_contig_types(types)
            )
            self.assertEqual(model.ploidy, 2)
            self.assertEqual(model.unitig_groups["C"], ("chr1_group1", "chr1_group2"))
            self.assertEqual(summary["validated"], 3)

    def test_group_selection_keeps_siblings_for_fair_assignment(self):
        model = self.make_model()
        args = SimpleNamespace(groups=["chr1_group1"], chromosomes=None)
        analysis, output = phase_workflow.select_group_models(model, args)
        self.assertEqual(set(analysis.groups), {"chr1_group1", "chr1_group2"})
        self.assertEqual(set(output.groups), {"chr1_group1"})
        self.assertEqual(
            analysis.unitig_groups["C"], ("chr1_group1", "chr1_group2")
        )

    def test_rerun_from_reuses_only_earlier_stage_checkpoints(self):
        args = SimpleNamespace(resume=True, rerun_from="assemble")
        self.assertTrue(phase_workflow.resume_stage(args, "assign"))
        self.assertTrue(phase_workflow.resume_stage(args, "extract"))
        self.assertFalse(phase_workflow.resume_stage(args, "assemble"))
        self.assertFalse(phase_workflow.resume_stage(args, "scaffold"))

        previous = {
            "stages": {
                "assemble": {
                    "parameters": {"hifiasm": "hifiasm"},
                    "runtime": {"jobs": 4, "threads_per_job": 10},
                }
            }
        }
        current = {
            "stages": {
                "assemble": {
                    "parameters": {"hifiasm": "hifiasm"},
                    "runtime": {"jobs": 1, "threads_per_job": 40},
                }
            }
        }
        phase_workflow.validate_resume_configuration(
            previous, current, ["assemble"]
        )

    def test_progress_logger_formats_large_counts(self):
        progress = phase_workflow.ProgressLog("test", every=1)
        with self.assertLogs("phap.phase_reads", level="INFO") as captured:
            progress(100_000, None, time.monotonic() - 1)
            progress(200_000, 200_000, time.monotonic() - 2)
        self.assertIn("100,000 processed", captured.output[0])
        self.assertIn("200,000/200,000", captured.output[1])

    def test_collapsed_reads_are_disjoint_conserved_and_depth_balanced(self):
        model = self.make_model()
        observations = {}
        for index, length in enumerate((1200, 1100, 1000, 900, 800, 700), 1):
            read_id = f"collapsed{index}"
            observations[read_id] = phase.ReadObservation(
                read_id,
                length,
                [self.hit("C", model.unitig_groups["C"], length=min(length, 1000))],
            )
        observations["unique"] = phase.ReadObservation(
            "unique", 1000, [self.hit("A", model.unitig_groups["A"])]
        )
        parameters = phase.AssignmentParameters(
            min_group_margin=0.10, min_total_aligned_bp=500
        )
        target_lengths = {"chr1_group1": 2000, "chr1_group2": 2000}
        decisions, summary = phase.assign_long_reads(
            observations, model, target_lengths, parameters, seed=100
        )
        self.assertEqual(set(decisions), set(observations))
        self.assertTrue(all(value.status == "assigned" for value in decisions.values()))
        self.assertEqual(decisions["unique"].group, "chr1_group1")
        collapsed = [
            value for key, value in decisions.items() if key.startswith("collapsed")
        ]
        self.assertTrue(
            all(value.basis == "collapsed_depth_balanced" for value in collapsed)
        )
        self.assertEqual(len({value.read_id for value in collapsed}), 6)
        group_bp = [
            summary["groups"][group]["read_bp"]
            for group in ("chr1_group1", "chr1_group2")
        ]
        self.assertLessEqual(abs(group_bp[0] - group_bp[1]), 1200)

        reversed_observations = dict(reversed(list(observations.items())))
        repeated, _ = phase.assign_long_reads(
            reversed_observations, model, target_lengths, parameters, seed=100
        )
        self.assertEqual(
            {key: value.group for key, value in decisions.items()},
            {key: value.group for key, value in repeated.items()},
        )

    def test_hic_pair_uses_both_mates_and_preserves_conflicts(self):
        model = self.make_model()
        pairs = {
            "compatible": {
                1: phase.ReadObservation(
                    "compatible", 150, [self.hit("A", model.unitig_groups["A"])]
                ),
                2: phase.ReadObservation(
                    "compatible", 150, [self.hit("C", model.unitig_groups["C"])]
                ),
            },
            "conflict": {
                1: phase.ReadObservation(
                    "conflict", 150, [self.hit("A", model.unitig_groups["A"])]
                ),
                2: phase.ReadObservation(
                    "conflict", 150, [self.hit("B", model.unitig_groups["B"])]
                ),
            },
            "collapsed": {
                1: phase.ReadObservation(
                    "collapsed", 150, [self.hit("C", model.unitig_groups["C"])]
                ),
                2: phase.ReadObservation(
                    "collapsed", 150, [self.hit("C", model.unitig_groups["C"])]
                ),
            },
        }
        decisions, _ = phase.assign_hic_pairs(
            pairs,
            model,
            {"chr1_group1": 2000, "chr1_group2": 2000},
            phase.AssignmentParameters(),
            seed=102,
        )
        self.assertEqual(decisions["compatible"].group, "chr1_group1")
        self.assertEqual(decisions["compatible"].basis, "hic_anchor_compatible_mate")
        self.assertEqual(decisions["conflict"].status, "unassigned")
        self.assertEqual(decisions["conflict"].basis, "hic_mates_conflict")
        self.assertEqual(decisions["collapsed"].status, "assigned")
        self.assertEqual(decisions["collapsed"].basis, "hic_collapsed_pair_balanced")

    def test_fastq_is_scanned_once_and_dispatched_exclusively(self):
        model = self.make_model()
        decisions = {
            "r1": phase.AssignmentDecision(
                "r1", 4, "assigned", "chr1_group1", "test", (), 1, 0, 1, 4, 1, ("A",)
            ),
            "r2": phase.AssignmentDecision(
                "r2", 4, "assigned", "chr1_group2", "test", (), 1, 0, 1, 4, 1, ("B",)
            ),
        }
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            fastq = work / "reads.fq"
            fastq.write_text(
                "@r1\nAAAA\n+\nIIII\n"
                "@r2\nCCCC\n+\nIIII\n"
                "@unused\nGGGG\n+\nIIII\n"
            )
            output = work / "out"
            summary = phase.dispatch_single_fastq(
                fastq, decisions, model.groups, output, "HiFi"
            )
            self.assertEqual(summary["full_input_scans"], 1)
            self.assertEqual(summary["written_reads"], 2)
            with gzip.open(output / "chr1_group1.HiFi.fq.gz", "rt") as handle:
                self.assertIn("@r1", handle.read())
            with gzip.open(output / "chr1_group2.HiFi.fq.gz", "rt") as handle:
                text = handle.read()
                self.assertIn("@r2", text)
                self.assertNotIn("@r1", text)

    def test_paired_fastq_normalizes_mate_suffixes(self):
        model = self.make_model()
        decisions = {
            "pair": phase.AssignmentDecision(
                "pair", 8, "assigned", "chr1_group1", "test", (), 1, 0, 1, 4, 2, ("A",)
            )
        }
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            read1 = work / "r1.fq"
            read2 = work / "r2.fq"
            read1.write_text("@pair/1 extra\nAAAA\n+\nIIII\n")
            read2.write_text("@pair/2 extra\nTTTT\n+\nIIII\n")
            summary = phase.dispatch_paired_fastq(
                read1, read2, decisions, model.groups, work / "out"
            )
            self.assertEqual(summary["written_pairs"], 1)
            self.assertEqual(summary["missing_pairs"], 0)

    @unittest.skipUnless(
        all(shutil.which(tool) for tool in ("seqkit", "pigz", "gawk")),
        "seqkit, pigz, and gawk are required for the fast backend test",
    )
    def test_seqkit_backend_preserves_reads_and_filters_quality(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            database = work / "assignments.sqlite"
            with sqlite3.connect(database) as connection:
                connection.execute(
                    "CREATE TABLE decisions ("
                    "read_id TEXT PRIMARY KEY, status TEXT, group_name TEXT)"
                )
                connection.executemany(
                    "INSERT INTO decisions VALUES (?, 'assigned', ?)",
                    (("keep", "chr1_group1"), ("filter", "chr1_group2")),
                )
                connection.execute(
                    "CREATE INDEX decisions_status_group "
                    "ON decisions(status, group_name)"
                )
            fastq = work / "reads.fq"
            fastq.write_text(
                "@keep description\nACGT\n+\n!I!I\n"
                "@filter\nTGCA\n+\n!!!!\n"
                "@unused\nGGGG\n+\nIIII\n"
            )
            output = work / "output"
            summary = phase.dispatch_single_fastq_seqkit(
                fastq, database,
                ("chr1_group1", "chr1_group2"), output, "ONT",
                seqkit=shutil.which("seqkit"),
                pigz=shutil.which("pigz"),
                gawk=shutil.which("gawk"),
                threads=2, min_length=1, min_mean_quality=20,
                progress_every=1,
            )
            self.assertEqual(summary["backend"], "seqkit")
            self.assertEqual(summary["input_records"], 3)
            self.assertEqual(summary["written_reads"], 1)
            self.assertEqual(summary["filtered_assigned_reads"], 1)
            self.assertEqual(summary["missing_reads"], 0)
            with gzip.open(output / "chr1_group1.ONT.fq.gz", "rt") as handle:
                self.assertEqual(
                    handle.read(), "@keep description\nACGT\n+\n!I!I\n"
                )
            with gzip.open(output / "chr1_group2.ONT.fq.gz", "rt") as handle:
                self.assertEqual(handle.read(), "")

    @unittest.skipUnless(
        all(shutil.which(tool) for tool in ("seqkit", "pigz", "gawk")),
        "seqkit, pigz, and gawk are required for fast-v1 Hi-C extraction",
    )
    def test_fast_v1_hic_seqkit_extraction_preserves_pairs(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            assignment = work / "assignment"
            assignment.mkdir()
            (assignment / "chr1_group1.read_ids.txt").write_text("pair1\npair3\n")
            (assignment / "chr1_group2.read_ids.txt").write_text("pair2\n")
            (assignment / "summary.json").write_text(json.dumps({
                "assignments": {
                    "groups": {
                        "chr1_group1": {"pairs": 2},
                        "chr1_group2": {"pairs": 1},
                    }
                }
            }))
            read1 = work / "hic1.fq"
            read2 = work / "hic2.fq"
            self.write_fastq(read1, ["pair1", "pair2", "pair3", "unused"], 20, 1)
            self.write_fastq(read2, ["pair1", "pair2", "pair3", "unused"], 20, 2)
            output = work / "output"
            summary = phase.dispatch_paired_fastq_fast_v1(
                read1, read2, assignment,
                ("chr1_group1", "chr1_group2"), output,
                seqkit=shutil.which("seqkit"),
                pigz=shutil.which("pigz"),
                gawk=shutil.which("gawk"),
                threads=2, jobs=2,
            )
            with gzip.open(output / "chr1_group1.Hi-C.1.fq.gz", "rt") as handle:
                text1 = handle.read()
            with gzip.open(output / "chr1_group1.Hi-C.2.fq.gz", "rt") as handle:
                text2 = handle.read()

        self.assertEqual(summary["written_pairs"], 3)
        self.assertEqual(summary["missing_pairs"], 0)
        self.assertIn("@pair1/1", text1)
        self.assertIn("@pair3/1", text1)
        self.assertNotIn("@pair2/1", text1)
        self.assertIn("@pair1/2", text2)
        self.assertIn("@pair3/2", text2)

    def write_bam(self, path, records, read_length):
        header = {
            "HD": {"VN": "1.6"},
            "SQ": [
                {"SN": "A", "LN": 2000},
                {"SN": "B", "LN": 2000},
                {"SN": "C", "LN": 2000},
            ],
        }
        with pysam.AlignmentFile(path, "wb", header=header) as bam:
            for record in records:
                name, reference, flag = record[:3]
                mapq = record[3] if len(record) > 3 else 60
                alignment = pysam.AlignedSegment()
                alignment.query_name = name
                alignment.query_sequence = "A" * read_length
                alignment.flag = flag
                alignment.reference_id = {"A": 0, "B": 1, "C": 2}[reference]
                alignment.reference_start = 0
                alignment.mapping_quality = mapq
                alignment.cigar = [(0, read_length)]
                mate_reference = record[4] if len(record) > 4 else reference
                alignment.next_reference_id = {
                    "A": 0, "B": 1, "C": 2
                }[mate_reference]
                alignment.next_reference_start = 0
                alignment.template_length = read_length * 2
                alignment.query_qualities = pysam.qualitystring_to_array("I" * read_length)
                alignment.set_tag("NM", 0)
                bam.write(alignment)

    def test_fast_v1_hic_is_disjoint_uses_rnext_and_ignores_mapq(self):
        model = self.make_model()
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            bam = work / "hic.bam"
            records = []
            for name, reference, mate_reference, mapq in (
                ("unique1", "A", "A", 60),
                ("unique2", "B", "B", 60),
                ("collapsed", "C", "C", 60),
                ("anchor", "A", "C", 60),
                ("conflict", "A", "B", 60),
                ("mapq0", "B", "B", 0),
            ):
                records.extend(
                    [
                        (name, reference, 99, mapq, mate_reference),
                        (name, mate_reference, 147, mapq, reference),
                    ]
                )
            self.write_bam(bam, records, 100)
            output = work / "fast"
            summary = phase.assign_hic_fast_v1(
                bam,
                output,
                model,
                phase.AssignmentParameters(
                    min_alignment_length=50,
                    min_identity=0.90,
                    collapsed_policy="balanced",
                ),
                seed=100,
            )
            group1 = set(
                (output / "chr1_group1.read_ids.txt").read_text().splitlines()
            )
            group2 = set(
                (output / "chr1_group2.read_ids.txt").read_text().splitlines()
            )

        self.assertFalse(group1 & group2)
        self.assertEqual(group1 | group2, {
            "unique1", "unique2", "collapsed", "anchor", "mapq0"
        })
        self.assertIn("anchor", group1)
        self.assertIn("mapq0", group2)
        self.assertNotIn("conflict", group1 | group2)
        self.assertEqual(summary["assignments"]["assigned"], 5)
        self.assertEqual(summary["assignments"]["unassigned"], 1)
        self.assertEqual(
            summary["assignments"]["basis_counts"]["hic_mates_conflict"], 1
        )
        self.assertFalse(any(output.glob("*.sqlite*")))

    @unittest.skipUnless(
        all(shutil.which(tool) for tool in ("seqkit", "pigz", "gawk")),
        "seqkit, pigz, and gawk are required for fast-v1 workflow",
    )
    def test_public_fast_v1_hic_workflow_uses_no_sqlite(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            groups = work / "groups.txt"
            groups.write_text("chr1_group1\tA C\nchr1_group2\tB C\n")
            types = work / "types.tsv"
            types.write_text(
                "contig_ID\taverage_depth\tcontig_type\n"
                "A\t20\thaplotig\nB\t20\thaplotig\nC\t40\tdiplotig\n"
            )
            bam = work / "hic.bam"
            records = []
            for name, reference, mate_reference in (
                ("p1", "A", "A"),
                ("p2", "B", "B"),
                ("pc", "C", "C"),
            ):
                records.extend([
                    (name, reference, 99, 0, mate_reference),
                    (name, mate_reference, 147, 0, reference),
                ])
            self.write_bam(bam, records, 100)
            hic1 = work / "hic1.fq"
            hic2 = work / "hic2.fq"
            self.write_fastq(hic1, ["p1", "p2", "pc"], 100, 1)
            self.write_fastq(hic2, ["p1", "p2", "pc"], 100, 2)
            output = work / "phase"
            completed = subprocess.run(
                [
                    sys.executable, str(PHASE_SCRIPT),
                    "--bam-hic", str(bam),
                    "--contig-type", str(types),
                    "--group", str(groups),
                    "--hic1", str(hic1),
                    "--hic2", str(hic2),
                    "--output-dir", str(output),
                    "--data-types", "hic",
                    "--stop-after", "extract",
                    "--hic-assignment-backend", "fast-v1",
                    "--extract-threads", "2",
                    "--jobs", "2",
                ],
                check=True, capture_output=True, text=True,
            )
            with (output / "02.reads" / "extraction_summary.json").open() as handle:
                summary = json.load(handle)
            sqlite_outputs = list(output.rglob("*.sqlite*"))

        self.assertIn("phase_reads completed", completed.stdout)
        self.assertEqual(summary["hic"]["written_pairs"], 3)
        self.assertEqual(summary["hic"]["missing_pairs"], 0)
        self.assertEqual(sqlite_outputs, [])

    def test_mapq_zero_primary_alignment_is_retained_and_scored(self):
        model = self.make_model()
        with tempfile.TemporaryDirectory() as temporary:
            bam = Path(temporary) / "mapq0.bam"
            self.write_bam(bam, [("low_mapq", "A", 0, 0)], 1000)
            observations, reference_lengths, stats = phase.collect_long_read_evidence(
                bam, model, phase.AssignmentParameters()
            )

        self.assertIn("low_mapq", observations)
        self.assertEqual(stats["retained_alignment_records"], 1)
        hit = observations["low_mapq"].hits[0]
        self.assertEqual(hit.mapq, 0)
        self.assertEqual(hit.weight, 1.0)
        evaluation = phase.evaluate_observation(
            observations["low_mapq"], model, phase.AssignmentParameters()
        )
        self.assertEqual(evaluation.confident_group, "chr1_group1")
        self.assertEqual(evaluation.reason, "alignment_confident")

    def write_fastq(self, path, names, read_length, mate=None):
        with path.open("w") as handle:
            for name in names:
                suffix = f"/{mate}" if mate else ""
                handle.write(
                    f"@{name}{suffix}\n"
                    + "A" * read_length
                    + "\n+\n"
                    + "I" * read_length
                    + "\n"
                )

    def test_public_workflow_runs_assignment_and_one_pass_extraction(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            groups = work / "groups.txt"
            groups.write_text("chr1_group1\tA C\nchr1_group2\tB C\n")
            types = work / "types.tsv"
            types.write_text(
                "contig_ID\taverage_depth\tcontig_type\n"
                "A\t20\thaplotig\nB\t20\thaplotig\nC\t40\tdiplotig\n"
            )
            hifi_names = ["h1", "h2", "hc1", "hc2"]
            ont_names = ["o1", "o2", "oc1", "oc2"]
            hifi_bam = work / "hifi.bam"
            ont_bam = work / "ont.bam"
            hic_bam = work / "hic.bam"
            self.write_bam(
                hifi_bam,
                [("h1", "A", 0), ("h2", "B", 0), ("hc1", "C", 0), ("hc2", "C", 0)],
                1000,
            )
            self.write_bam(
                ont_bam,
                [("o1", "A", 0), ("o2", "B", 0), ("oc1", "C", 0), ("oc2", "C", 0)],
                1000,
            )
            hic_records = []
            for name, reference in (("p1", "A"), ("p2", "B"), ("pc", "C")):
                hic_records.extend([(name, reference, 99), (name, reference, 147)])
            self.write_bam(hic_bam, hic_records, 100)
            hifi = work / "hifi.fq"
            ont = work / "ont.fq"
            hic1 = work / "hic1.fq"
            hic2 = work / "hic2.fq"
            self.write_fastq(hifi, hifi_names, 1000)
            self.write_fastq(ont, ont_names, 1000)
            self.write_fastq(hic1, ["p1", "p2", "pc"], 100, mate=1)
            self.write_fastq(hic2, ["p1", "p2", "pc"], 100, mate=2)
            output = work / "phase"
            temp_root = work / "phase_tmp"
            base_command = [
                sys.executable,
                str(PHASE_SCRIPT),
                "--bam-hifi", str(hifi_bam),
                "--bam-ont", str(ont_bam),
                "--bam-hic", str(hic_bam),
                "--contig-type", str(types),
                "--group", str(groups),
                "--output-dir", str(output),
                "--temp-dir", str(temp_root),
                "--hic-assignment-backend", "sqlite",
            ]
            input_signatures = {
                path: (path.stat().st_size, path.stat().st_mtime_ns)
                for path in (hifi_bam, ont_bam, hic_bam, hifi, ont, hic1, hic2)
            }
            subprocess.run(
                base_command + ["--stop-after", "assign"],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                base_command
                + [
                    "--hifi", str(hifi),
                    "--ont", str(ont),
                    "--hic1", str(hic1),
                    "--hic2", str(hic2),
                    "--resume",
                    "--stop-after", "extract",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            with (output / "run_manifest.json").open() as handle:
                manifest = json.load(handle)
            self.assertEqual(manifest["completed_stages"], ["assign", "extract"])
            with (output / "02.reads" / "extraction_summary.json").open() as handle:
                summary = json.load(handle)
            self.assertEqual(summary["hifi"]["written_reads"], 4)
            self.assertEqual(summary["ont"]["written_reads"], 4)
            self.assertEqual(summary["hic"]["written_pairs"], 3)
            self.assertEqual(summary["hifi"]["full_input_scans"], 1)
            self.assertIn("[stage 1/4] assign", (output / "phase_reads.log").read_text())
            self.assertEqual(
                manifest["configuration"]["runtime"]["temp_dir"], str(temp_root)
            )
            self.assertEqual(list(temp_root.iterdir()), [])
            self.assertEqual(
                input_signatures,
                {
                    path: (path.stat().st_size, path.stat().st_mtime_ns)
                    for path in input_signatures
                },
            )

            hifi_only = work / "hifi_only"
            subprocess.run(
                [
                    sys.executable, str(PHASE_SCRIPT),
                    "--bam-hifi", str(hifi_bam),
                    "--contig-type", str(types),
                    "--group", str(groups),
                    "--output-dir", str(hifi_only),
                    "--data-types", "hifi",
                    "--stop-after", "assign",
                ],
                check=True, capture_output=True, text=True,
            )
            with (hifi_only / "01.assignments" / "assignment_summary.json").open() as handle:
                hifi_summary = json.load(handle)
            self.assertEqual(hifi_summary["data_types"], ["hifi"])
            self.assertTrue(
                (hifi_only / "01.assignments" / "hifi.assignments.sqlite").exists()
            )
            self.assertFalse(
                (hifi_only / "01.assignments" / "ont.assignments.sqlite").exists()
            )

    def test_public_workflow_no_collapse_needs_no_contig_type(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            groups = work / "groups.txt"
            groups.write_text("chr1_group1\tA\nchr1_group2\tB\n")
            hifi_bam = work / "hifi.bam"
            self.write_bam(
                hifi_bam,
                [("h1", "A", 0), ("h2", "B", 0)],
                1000,
            )
            output = work / "phase"
            subprocess.run(
                [
                    sys.executable, str(PHASE_SCRIPT),
                    "--bam-hifi", str(hifi_bam),
                    "--group", str(groups),
                    "--output-dir", str(output),
                    "--data-types", "hifi",
                    "--stop-after", "assign",
                    "--no-collapse",
                ],
                check=True, capture_output=True, text=True,
            )
            summary = json.loads(
                (output / "01.assignments" / "assignment_summary.json").read_text()
            )
            validation = summary["group_model"]["dosage_validation"]
            self.assertEqual(validation["mode"], "single_copy")
            self.assertEqual(validation["validated"], 2)

    def test_disk_evidence_checkpoint_resumes_from_bam_offset(self):
        model = self.make_model()
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            bam = work / "reads.bam"
            self.write_bam(
                bam,
                [("r1", "A", 0), ("r2", "B", 0), ("r3", "C", 0)],
                1000,
            )
            database = work / "evidence.sqlite"
            calls = 0

            def interrupt(processed, total, started, detail=None):
                nonlocal calls
                calls += 1
                if calls == 1:
                    raise RuntimeError("simulated interruption")

            with self.assertRaises(RuntimeError):
                phase.collect_evidence_disk(
                    bam, model, phase.AssignmentParameters(), database,
                    "checkpoint", progress=interrupt, progress_every=1,
                )
            result = phase.collect_evidence_disk(
                bam, model, phase.AssignmentParameters(), database,
                "checkpoint", progress_every=1,
            )
            self.assertGreaterEqual(result["resumed_from_record"], 1)
            with sqlite3.connect(database) as connection:
                retained = connection.execute(
                    "SELECT COUNT(*) FROM evidence"
                ).fetchone()[0]
            self.assertEqual(retained, 3)

    def test_assembly_stage_uses_bounded_jobs_and_validates_output(self):
        model = self.make_model()
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            reads = work / "reads"
            reads.mkdir()
            for group in model.groups:
                for suffix in ("HiFi", "ONT"):
                    with gzip.open(reads / f"{group}.{suffix}.fq.gz", "wt") as handle:
                        handle.write("@read\nAAAA\n+\nIIII\n")
            fake = work / "fake_hifiasm.py"
            calls = work / "hifiasm.calls"
            fake.write_text(
                "#!/usr/bin/env python3\n"
                "import pathlib, sys\n"
                "prefix = sys.argv[sys.argv.index('-o') + 1]\n"
                f"count = pathlib.Path({str(calls)!r})\n"
                "with count.open('a') as handle: handle.write(prefix + '\\n')\n"
                "pathlib.Path(prefix + '.bp.p_ctg.gfa').write_text('S\\tctg\\tAAAA\\n')\n"
            )
            fake.chmod(0o755)
            output = work / "assembly"
            args = SimpleNamespace(
                hifiasm=str(fake), threads_per_job=2, jobs=2
            )
            summary = phase_workflow.run_assembly_stage(
                args, model, reads, output
            )
            self.assertEqual(summary["groups"], 2)
            for group in model.groups:
                self.assertTrue(
                    (output / group / f"{group}.asm.bp.p_ctg.gfa").exists()
                )
            self.assertEqual(len(calls.read_text().splitlines()), 2)
            missing_group = "chr1_group2"
            (output / missing_group / f"{missing_group}.asm.bp.p_ctg.gfa").unlink()
            args.resume = True
            phase_workflow.run_assembly_stage(args, model, reads, output)
            self.assertEqual(len(calls.read_text().splitlines()), 3)

    def test_failed_assembly_preserves_previous_stage_output(self):
        model = self.make_model()
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            reads = work / "reads"
            reads.mkdir()
            for group in model.groups:
                with gzip.open(reads / f"{group}.HiFi.fq.gz", "wt") as handle:
                    handle.write("@read\nAAAA\n+\nIIII\n")
                with gzip.open(reads / f"{group}.ONT.fq.gz", "wt") as handle:
                    pass
            fake = work / "failed_hifiasm.py"
            fake.write_text("#!/usr/bin/env python3\nraise SystemExit(3)\n")
            fake.chmod(0o755)
            output = work / "assembly"
            output.mkdir()
            (output / "keep.txt").write_text("previous")
            args = SimpleNamespace(
                hifiasm=str(fake), threads_per_job=1, jobs=1
            )
            with self.assertRaises(subprocess.CalledProcessError):
                phase_workflow.run_assembly_stage(args, model, reads, output)
            self.assertEqual((output / "keep.txt").read_text(), "previous")

    def test_scaffold_helpers_convert_gfa_and_propagate_pipe_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            gfa = work / "assembly.gfa"
            fasta = work / "assembly.fa"
            gfa.write_text("H\tVN:Z:1.0\nS\tctg1\tACGT\nS\tctg2\tTTAA\n")
            records = phase_workflow.gfa_to_fasta(gfa, fasta)
            self.assertEqual(records, 2)
            self.assertEqual(fasta.read_text(), ">ctg1\nACGT\n>ctg2\nTTAA\n")

            commands = [
                [sys.executable, "-c", "print('alignment')"],
                [
                    sys.executable,
                    "-c",
                    "import sys; sys.stdin.buffer.read(); raise SystemExit(7)",
                ],
            ]
            with self.assertRaises(subprocess.CalledProcessError) as caught:
                phase_workflow.run_pipe(commands, work, work / "mapping_pipe")
            self.assertEqual(caught.exception.returncode, 7)

    def test_scaffold_only_uses_external_inputs_and_group_checkpoints(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            groups = work / "groups.txt"
            groups.write_text("chr1_group1\tA\nchr1_group2\tB\n")
            assemblies = work / "assemblies"
            assembly_group = assemblies / "chr1_group1"
            assembly_group.mkdir(parents=True)
            (assembly_group / "chr1_group1.asm.bp.p_ctg.gfa").write_text(
                "H\tVN:Z:1.0\nS\tctg1\tACGTACGT\n"
            )
            reads = work / "hic_reads"
            reads.mkdir()
            for mate in (1, 2):
                with gzip.open(
                    reads / f"chr1_group1.Hi-C.{mate}.fq.gz", "wt"
                ) as handle:
                    handle.write(f"@pair/{mate}\nACGT\n+\nIIII\n")

            tools = work / "tools"
            tools.mkdir()

            def executable(name, text):
                path = tools / name
                path.write_text("#!/usr/bin/env python3\n" + text)
                path.chmod(0o755)
                return path

            bwa = executable(
                "bwa",
                "import pathlib, sys\n"
                "if sys.argv[1] == 'index':\n"
                "    pathlib.Path(sys.argv[2] + '.bwt').write_text('index')\n"
                "elif sys.argv[1] == 'mem':\n"
                "    sys.stdout.write('@HD\\tVN:1.6\\n')\n",
            )
            samblaster = executable(
                "samblaster",
                "import sys\n"
                "sys.stdout.buffer.write(sys.stdin.buffer.read())\n",
            )
            samtools = executable(
                "samtools",
                "import pathlib, sys\n"
                "data = sys.stdin.buffer.read()\n"
                "output = pathlib.Path(sys.argv[sys.argv.index('-o') + 1])\n"
                "output.write_bytes(data or b'BAM')\n",
            )
            filter_bam = executable(
                "filter_bam",
                "import sys\n"
                "sys.stdout.buffer.write(b'FILTERED')\n",
            )
            haphic = executable(
                "haphic",
                "import pathlib, sys\n"
                "build = pathlib.Path('04.build')\n"
                "build.mkdir()\n"
                "(build / 'scaffolds.fa').write_text('>scaffold\\nACGT\\n')\n"
                "(build / 'juicebox.sh').write_text(\n"
                "    '#!/bin/bash\\n'\n"
                "    f'test -s {sys.argv[2]}\\n'\n"
                "    f'test -s {sys.argv[3]}\\n'\n"
                ")\n",
            )
            output = work / "scaffold_run"
            temp_root = work / "temp"
            command = [
                sys.executable, str(PHASE_SCRIPT),
                "--scaffold-only",
                "--group", str(groups),
                "--assembly-dir", str(assemblies),
                "--hic-reads-dir", str(reads),
                "--output-dir", str(output),
                "--temp-dir", str(temp_root),
                "--groups", "chr1_group1",
                "--jobs", "1",
                "--threads-per-job", "1",
                "--bwa", str(bwa),
                "--samtools", str(samtools),
                "--samblaster", str(samblaster),
                "--filter-bam", str(filter_bam),
                "--haphic", str(haphic),
            ]
            first = subprocess.run(
                command, check=True, capture_output=True, text=True
            )
            juicebox = (
                output / "04.scaffold" / "chr1_group1"
                / "02.haphic" / "04.build" / "juicebox.sh"
            )
            stale_group = temp_root / "old" / ".chr1_group1.temporary"
            juicebox.write_text(
                "#!/bin/bash\n"
                f"test -s {stale_group}/01.hic_mapping/chr1_group1.p_ctg.fa\n"
                f"test -s {stale_group}/01.hic_mapping/HiC.filtered.bam\n"
            )
            second = subprocess.run(
                command + ["--resume"],
                check=True, capture_output=True, text=True,
            )
            (assembly_group / "chr1_group1.asm.bp.p_ctg.gfa").write_text(
                "H\tVN:Z:1.0\nS\tctg1\tACGTACGTAA\n"
            )
            changed_input = subprocess.run(
                command + ["--resume"], capture_output=True, text=True
            )
            with (output / "run_manifest.json").open() as handle:
                manifest = json.load(handle)
            build = (
                output / "04.scaffold" / "chr1_group1"
                / "02.haphic" / "04.build" / "scaffolds.fa"
            )
            checkpoint = (
                output / "04.scaffold" / ".checkpoints"
                / "chr1_group1.complete.json"
            )
            build_exists = build.exists()
            checkpoint_exists = checkpoint.exists()
            juicebox_text = juicebox.read_text()
            transient_paths_remain = str(temp_root.resolve()) in juicebox_text
            juicebox_run = subprocess.run(
                ["bash", str(juicebox)], cwd=juicebox.parent,
                capture_output=True, text=True,
            )

        self.assertIn("completed stages: scaffold", first.stdout)
        self.assertIn("completed stages: scaffold", second.stdout)
        self.assertEqual(manifest["mode"], "scaffold-only")
        self.assertEqual(manifest["completed_stages"], ["scaffold"])
        self.assertEqual(
            set(manifest["configuration"]["stages"]["scaffold"]["inputs"]["groups"]),
            {"chr1_group1"},
        )
        self.assertTrue(build_exists)
        self.assertTrue(checkpoint_exists)
        self.assertFalse(transient_paths_remain)
        self.assertIn("../../01.hic_mapping/chr1_group1.p_ctg.fa", juicebox_text)
        self.assertIn("../../01.hic_mapping/HiC.filtered.bam", juicebox_text)
        self.assertEqual(juicebox_run.returncode, 0, juicebox_run.stderr)
        self.assertIn("repaired 2 stale juicebox path(s)", second.stderr)
        self.assertIn("checkpoint complete; skipped chr1_group1", second.stderr)
        self.assertNotEqual(changed_input.returncode, 0)
        self.assertIn(
            "resume parameters for completed stage scaffold differ",
            changed_input.stderr,
        )


if __name__ == "__main__":
    unittest.main()
