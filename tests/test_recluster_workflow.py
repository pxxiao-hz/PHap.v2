from __future__ import annotations

import os
import pickle
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ReclusterWorkflowTests(unittest.TestCase):
    def test_chromosome_recluster_uses_selected_score_mode_and_dynamic_groups(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fasta = directory / "scaffold.putg.fa.copy.putg.fa"
            fasta.write_text(
                ">anchor1\nGATC\n>anchor2\nGATC\n>candidate\nGATC\n",
                encoding="utf-8",
            )
            dosage = directory / "dosage.tsv"
            dosage.write_text(
                "contig_ID\taverage_depth\tcontig_type\tdosage\tstatus\n"
                "anchor1\t20\thaplotig\t1\tassigned\n"
                "anchor2\t20\thaplotig\t1\tassigned\n"
                "candidate\t20\thaplotig\t1\tassigned\n",
                encoding="utf-8",
            )
            groups = directory / "groups.tsv"
            groups.write_text(
                "group1\t1\tanchor1\n"
                "group2\t1\tanchor2\n",
                encoding="utf-8",
            )
            links = directory / "links.pkl"
            with links.open("wb") as output:
                pickle.dump(
                    {
                        ("candidate", "anchor1"): 20,
                        ("candidate", "anchor2"): 2,
                    },
                    output,
                )
            clm = directory / "links.clm"
            clm.write_text("", encoding="utf-8")
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "utils.chr_uncluster_recluster",
                    "--fasta",
                    str(fasta),
                    "--contig_type",
                    str(dosage),
                    "--full_links",
                    str(links),
                    "--clusters_file",
                    str(groups),
                    "--clm_file",
                    str(clm),
                    "--ploidy",
                    "2",
                    "--hic_score_mode",
                    "raw",
                    "--min_hic_score",
                    "1",
                    "--min_hic_margin",
                    "1",
                ],
                cwd=directory,
                env={**os.environ, "PYTHONPATH": str(ROOT)},
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            cluster_text = (directory / "group.reassignment.cluster.txt").read_text(
                encoding="utf-8"
            )
            audit_text = (directory / "recluster_scores.tsv").read_text(
                encoding="utf-8"
            )

        self.assertIn(
            "scaffold.putg.fa.copy_group1\tanchor1\tcandidate",
            cluster_text,
        )
        self.assertIn("candidate\tgroup1\t20.000000", audit_text)
        self.assertIn("candidate\tgroup2\t2.000000", audit_text)


if __name__ == "__main__":
    unittest.main()
