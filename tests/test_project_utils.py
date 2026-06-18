import tempfile
import unittest
from pathlib import Path

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from project_utils import cpg_count, gc_content, read_sequence_file, revcomp, validate_seq


class ProjectUtilsTests(unittest.TestCase):
    def test_sequence_validation_and_metrics(self):
        seq = "ACGT" * 50
        self.assertTrue(validate_seq(seq))
        self.assertEqual(revcomp(seq), seq)
        self.assertAlmostEqual(gc_content(seq), 0.5)
        self.assertEqual(cpg_count("ACG" * 66 + "AA"), 66)

    def test_read_sequence_file_plain_text(self):
        seq = "A" * 200
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "seqs.txt"
            path.write_text(seq + "\ninvalid\n")
            self.assertEqual(read_sequence_file(path), [seq])

    def test_read_sequence_file_tsv(self):
        seq = "C" * 200
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "seqs.tsv"
            path.write_text("id\tsequence\nx\t" + seq + "\n")
            self.assertEqual(read_sequence_file(path), [seq])


if __name__ == "__main__":
    unittest.main()
