#!/usr/bin/env python3
"""Create the public-data manifest for the HSC/T-NK promoter project."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from project_utils import ensure_dir, write_json


SEED_MANIFEST = [
    {
        "source": "FANTOM5",
        "assays": "CAGE/promoter atlas",
        "use": "TSS choice and transcription-initiation signal",
        "cell_states": "blood and immune cell samples where available",
        "assembly": "hg19/source-dependent; lift to hg38",
        "url": "https://fantom.gsc.riken.jp/5/",
        "access": "public",
        "priority": "required",
    },
    {
        "source": "ENCODE",
        "assays": "ATAC/DNase/RAMPAGE/CAGE/RNA-seq",
        "use": "accessibility, initiation, expression, cCRE overlap",
        "cell_states": "CD34/HSPC and mature immune cells where available",
        "assembly": "hg38 preferred",
        "url": "https://www.encodeproject.org/",
        "access": "public",
        "priority": "required",
    },
    {
        "source": "BLUEPRINT",
        "assays": "RNA-seq, DNase/ATAC, histone marks",
        "use": "hematopoietic lineage expression and epigenome profiles",
        "cell_states": "HSC/HSPC, progenitors, myeloid, lymphoid",
        "assembly": "source-dependent; harmonize to hg38",
        "url": "https://projects.ensembl.org/blueprint/",
        "access": "mixed public/controlled",
        "priority": "required-if-accessible",
    },
    {
        "source": "Roadmap Epigenomics",
        "assays": "DNase, histone marks, RNA-seq",
        "use": "fallback hematopoietic epigenome coverage",
        "cell_states": "CD34 and immune/blood reference states",
        "assembly": "hg19; lift to hg38",
        "url": "https://egg2.wustl.edu/roadmap/web_portal/processed_data.html",
        "access": "public",
        "priority": "high",
    },
    {
        "source": "DICE",
        "assays": "RNA-seq and immune epigenomics",
        "use": "mature immune off-target and T/NK reference profiles",
        "cell_states": "T, NK, B, monocytes and immune subsets",
        "assembly": "source-dependent; harmonize to hg38",
        "url": "https://dice-database.org/",
        "access": "public",
        "priority": "high",
    },
    {
        "source": "Human Cell Atlas / 10x",
        "assays": "scRNA/scATAC/multiome",
        "use": "CD34 differentiation pseudobulk profiles",
        "cell_states": "HSC/HSPC, progenitors, mature bone marrow lineages",
        "assembly": "source-dependent; harmonize to hg38/gene IDs",
        "url": "https://data.humancellatlas.org/",
        "access": "public",
        "priority": "high",
    },
    {
        "source": "ENCODE SCREEN",
        "assays": "cCRE annotations",
        "use": "promoter/enhancer annotation and blacklist-like filtering",
        "cell_states": "all available",
        "assembly": "hg38",
        "url": "https://screen.encodeproject.org/",
        "access": "public",
        "priority": "required",
    },
    {
        "source": "JASPAR/HOCOMOCO",
        "assays": "TF motif models",
        "use": "motif enrichment, unwanted motif filtering, interpretability",
        "cell_states": "not cell-state-specific",
        "assembly": "not applicable",
        "url": "https://jaspar.elixir.no/",
        "access": "public",
        "priority": "required",
    },
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default=Path(__file__).resolve().parents[1] / "reference_data" / "hematopoietic_manifest")
    args = parser.parse_args()
    out = ensure_dir(args.out_dir)
    csv_path = out / "public_hematopoietic_data_manifest.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(SEED_MANIFEST[0].keys()))
        writer.writeheader()
        writer.writerows(SEED_MANIFEST)
    write_json(out / "public_hematopoietic_data_manifest.json", SEED_MANIFEST)
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
