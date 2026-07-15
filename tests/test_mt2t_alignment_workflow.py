from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path
from typing import Any, Sequence
from unittest.mock import patch

from utils.mt2t_alignment_workflow import build_pairwise_overlap_paf


class Mt2tAlignmentWorkflowTests(unittest.TestCase):
    def test_two_contigs_use_argument_lists_and_safe_staging_names(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "input.fa"
            fasta.write_text(
                ">../unsafe;A description\nAAAACCCC\n>B\nCCCCGGGG\n",
                encoding="utf-8",
            )
            calls: list[tuple[str, ...]] = []

            def fake_run(
                command: Sequence[object],
                **kwargs: Any,
            ) -> subprocess.CompletedProcess[bytes]:
                args = tuple(str(value) for value in command)
                calls.append(args)
                if args[:2] == ("mash", "dist"):
                    return subprocess.CompletedProcess(
                        args,
                        0,
                        stdout=b"left\tright\t0.01\t0\t10/10\n",
                        stderr=b"",
                    )
                if args[0] == "minimap2":
                    kwargs["stdout"].write(
                        "B\t8\t0\t4\t+\t../unsafe;A\t8\t4\t8"
                        "\t4\t4\t60\ttp:A:P\n"
                    )
                return subprocess.CompletedProcess(args, 0, stdout=b"", stderr=b"")

            with patch("utils.mt2t_alignment_workflow.run_command", fake_run):
                merge = build_pairwise_overlap_paf(
                    fasta_path=str(fasta),
                    output_directory=str(root / "out"),
                    min_contig_length=1,
                    min_distance=0.15,
                    threads_per_alignment=1,
                    max_processes=1,
                    cpu_budget=1,
                    tool_versions=(("mash", "test"), ("minimap2", "test")),
                )
            merge_text = merge.read_text(encoding="utf-8")
            id_map = (root / "out" / "mt2t_split_id_map.tsv").read_text(
                encoding="utf-8"
            )
            manifest = (root / "out" / "mt2t_upstream_manifest.tsv").read_text(
                encoding="utf-8"
            )
        self.assertIn("../unsafe;A", merge_text)
        self.assertIn("contig_000001\t../unsafe;A", id_map)
        command_text = "\n".join(" ".join(call) for call in calls)
        self.assertNotIn("../unsafe;A", command_text)
        self.assertTrue(any(call[0] == "minimap2" for call in calls))
        self.assertIn("$STAGING", manifest)
        self.assertNotIn(".mt2t-alignment-", manifest)


if __name__ == "__main__":
    unittest.main()
