import math
import tempfile
import unittest
from pathlib import Path

import pandas as pd

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_real_stage2_inputs import (  # noqa: E402
    choose_tss_sources,
    expression_activity,
    normalize_activity_value,
    stable_promoter_id,
    validate_activity_rows,
)
from filter_and_rank_candidates import max_homopolymer  # noqa: E402


class RealStage2InputTests(unittest.TestCase):
    def test_stable_promoter_id_is_deterministic(self):
        a = stable_promoter_id("chr1", 10, 210, "+", "GENE1")
        b = stable_promoter_id("chr1", 10, 210, "+", "GENE1")
        c = stable_promoter_id("chr1", 11, 211, "+", "GENE1")
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)
        self.assertTrue(a.startswith("prom_"))

    def test_fantom_preferred_over_gencode_duplicate_gene(self):
        fantom = pd.DataFrame(
            [
                {
                    "chr": "chr1",
                    "tss0": 100,
                    "strand": "+",
                    "gene_id": "G1",
                    "gene_name": "GENE1",
                    "source": "FANTOM5",
                    "source_priority": 1,
                }
            ]
        )
        gencode = pd.DataFrame(
            [
                {
                    "chr": "chr1",
                    "tss0": 110,
                    "strand": "+",
                    "gene_id": "G1",
                    "gene_name": "GENE1",
                    "source": "GENCODE",
                    "source_priority": 2,
                }
            ]
        )
        chosen = choose_tss_sources([gencode, fantom])
        self.assertEqual(len(chosen), 1)
        self.assertEqual(chosen.iloc[0]["source"], "FANTOM5")
        self.assertEqual(int(chosen.iloc[0]["tss0"]), 100)

    def test_activity_validation_harmonizes_and_log_transforms(self):
        df = pd.DataFrame(
            {
                "promoter_id": ["p1"],
                "cell_type": ["CD4_T"],
                "assay": ["RNA-seq"],
                "value": [9.0],
            }
        )
        out = validate_activity_rows(df)
        self.assertEqual(out.iloc[0]["cell_type"], "T")
        self.assertEqual(out.iloc[0]["assay"], "expression")
        self.assertAlmostEqual(out.iloc[0]["value"], math.log1p(9.0))

    def test_expression_activity_joins_gene_ids(self):
        with tempfile.TemporaryDirectory() as td:
            expr = Path(td) / "expr.tsv"
            pd.DataFrame(
                {
                    "cell_type": ["NK", "HSC"],
                    "gene_id": ["G1", "G2"],
                    "value": [10.0, 5.0],
                    "source": ["unit", "unit"],
                }
            ).to_csv(expr, sep="\t", index=False)
            promoters = pd.DataFrame(
                {
                    "promoter_id": ["p1", "p2"],
                    "gene_id": ["G1", "G3"],
                    "gene_name": ["GENE1", "GENE3"],
                }
            )
            activity, _sources = expression_activity(promoters, expr)
            self.assertEqual(len(activity), 1)
            self.assertEqual(activity.iloc[0]["promoter_id"], "p1")
            self.assertEqual(activity.iloc[0]["cell_type"], "NK")

    def test_homopolymer_counter(self):
        self.assertEqual(max_homopolymer("AAACCCCCGT"), 5)
        self.assertEqual(max_homopolymer("ACGT"), 1)


if __name__ == "__main__":
    unittest.main()
