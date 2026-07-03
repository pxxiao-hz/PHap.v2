from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from typing import Any

try:
    import pysam  # type: ignore[import-not-found]
except ImportError:  # pragma: no cover - exercised in dependency-minimal local checks
    pysam = None  # type: ignore[assignment]


@unittest.skipUnless(pysam is not None, "pysam is required for BAM integration tests")
class BamReadEvidenceTests(unittest.TestCase):
    def test_hic_mates_are_combined_and_filtered_records_remain_audited(self) -> None:
        from utils.phase_reads_assemble_anchor import parse_bam_read_unitigs

        assert pysam is not None
        with tempfile.TemporaryDirectory() as temporary_directory:
            bam_path = Path(temporary_directory) / "reads.bam"
            header = {
                "HD": {"VN": "1.6", "SO": "unsorted"},
                "SQ": [
                    {"SN": "utgA", "LN": 100},
                    {"SN": "utgB", "LN": 100},
                ],
            }
            with pysam.AlignmentFile(bam_path, "wb", header=header) as output:
                output.write(self._alignment("pair1/1", reference_id=0, flag=65))
                output.write(self._alignment("pair1/2", reference_id=1, flag=129))
                output.write(self._alignment("qcfail", reference_id=0, flag=512))
                output.write(self._alignment("secondary", reference_id=0, flag=256))

            evidence = parse_bam_read_unitigs(
                str(bam_path),
                min_mapq=1,
                paired=True,
            )

        self.assertEqual(evidence["pair1"], {"utgA", "utgB"})
        self.assertEqual(evidence["qcfail"], set())
        self.assertEqual(evidence["secondary"], set())

    @staticmethod
    def _alignment(query_name: str, *, reference_id: int, flag: int) -> Any:
        assert pysam is not None
        alignment = pysam.AlignedSegment()
        alignment.query_name = query_name
        alignment.query_sequence = "A" * 20
        alignment.flag = flag
        alignment.reference_id = reference_id
        alignment.reference_start = 0
        alignment.mapping_quality = 60
        alignment.cigartuples = [(0, 20)]
        alignment.query_qualities = pysam.qualitystring_to_array("I" * 20)
        return alignment


if __name__ == "__main__":
    unittest.main()
