import sys
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from download_hca_stage2_expression import choose_h5ad_asset  # noqa: E402
from download_hpa_stage2_expression import hpa_to_expression_long  # noqa: E402


class PublicExpressionDownloaderTests(unittest.TestCase):
    def test_hca_asset_selection_prefers_public_h5ad(self):
        collection = {
            "datasets": [
                {
                    "dataset_id": "d1",
                    "title": "unrelated",
                    "cell_count": 100,
                    "assets": [{"filetype": "H5AD", "url": "https://example.org/a.h5ad", "filesize": 1}],
                },
                {
                    "dataset_id": "d2",
                    "title": "A Balanced Bone Marrow Reference Map of Hematopoietic Development",
                    "cell_count": 250000,
                    "disease": [{"label": "normal"}],
                    "assets": [{"filetype": "H5AD", "url": "https://example.org/b.h5ad", "filesize": 2}],
                },
            ]
        }
        dataset, asset = choose_h5ad_asset(collection)
        self.assertEqual(dataset["dataset_id"], "d2")
        self.assertEqual(asset["url"], "https://example.org/b.h5ad")

    def test_hpa_mapping_covers_required_lineages(self):
        df = pd.DataFrame(
            {
                "Gene": ["ENSG1", "ENSG1", "ENSG1", "ENSG1"],
                "Gene name": ["G1", "G1", "G1", "G1"],
                "Cell type": ["hematopoietic stem cells", "t-cells", "nk-cells", "megakaryocyte progenitors"],
                "nCPM": [1.0, 2.0, 3.0, 4.0],
            }
        )
        out, summary = hpa_to_expression_long(df, "unit", "unit.zip")
        self.assertTrue({"HSC", "HSPC", "T", "NK", "MEGAKARYOCYTE"} <= set(out["cell_type"]))
        self.assertEqual(summary["genes"].sum(), len(out))


if __name__ == "__main__":
    unittest.main()
