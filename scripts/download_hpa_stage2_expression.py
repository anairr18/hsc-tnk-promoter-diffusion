#!/usr/bin/env python3
"""Build expression_long.tsv from Human Protein Atlas single-cell type RNA."""

from __future__ import annotations

import argparse
import io
import zipfile
from pathlib import Path

import pandas as pd
import requests


HPA_SINGLE_CELL_URL = "https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip"

CELL_TYPE_MAP = {
    "hematopoietic stem cells": ["HSC", "HSPC"],
    "t-cells": ["T"],
    "thymocytes": ["T"],
    "nk-cells": ["NK"],
    "b-cells": ["B"],
    "monocytes": ["MYELOID"],
    "monocyte progenitors": ["MYELOID"],
    "macrophages": ["MYELOID"],
    "neutrophils": ["MYELOID"],
    "neutrophil progenitors": ["MYELOID"],
    "erythrocytes": ["ERYTHROID"],
    "erythrocyte progenitors": ["ERYTHROID"],
    "megakaryocyte-erythroid progenitors": ["ERYTHROID", "MEGAKARYOCYTE"],
    "megakaryocytes": ["MEGAKARYOCYTE"],
    "megakaryocyte progenitors": ["MEGAKARYOCYTE"],
}


def download_bytes(url: str) -> bytes:
    response = requests.get(url, timeout=180)
    response.raise_for_status()
    return response.content


def hpa_to_expression_long(df: pd.DataFrame, source: str, accession: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"Gene", "Gene name", "Cell type", "nCPM"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"HPA table missing expected columns: {sorted(missing)}")
    rows = []
    summary_rows = []
    for hpa_cell, harmonized_cells in CELL_TYPE_MAP.items():
        sub = df[df["Cell type"].astype(str).str.lower() == hpa_cell]
        if sub.empty:
            continue
        for cell_type in harmonized_cells:
            summary_rows.append({"cell_type": cell_type, "hpa_cell_type": hpa_cell, "genes": int(len(sub))})
            for _, row in sub.iterrows():
                rows.append(
                    {
                        "cell_type": cell_type,
                        "gene_id": str(row["Gene"]).split(".")[0],
                        "gene_name": row["Gene name"],
                        "value": float(row["nCPM"]),
                        "source": source,
                        "accession": accession,
                        "replicate": hpa_cell,
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(summary_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default=HPA_SINGLE_CELL_URL)
    parser.add_argument("--download-dir", default=Path(__file__).resolve().parents[1] / "downloads" / "hpa")
    parser.add_argument("--output", default=Path(__file__).resolve().parents[1] / "data" / "hsc_tnk_real" / "hpa_expression_long.tsv")
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--source-label", default="HPA_SINGLE_CELL_TYPE")
    args = parser.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    dl_dir = Path(args.download_dir)
    dl_dir.mkdir(parents=True, exist_ok=True)
    zip_path = dl_dir / "rna_single_cell_type.tsv.zip"
    if args.force_download or not zip_path.exists():
        zip_path.write_bytes(download_bytes(args.url))

    with zipfile.ZipFile(zip_path) as zf:
        names = [n for n in zf.namelist() if n.endswith(".tsv")]
        if not names:
            raise SystemExit("No TSV found in HPA zip.")
        with zf.open(names[0]) as handle:
            df = pd.read_csv(io.TextIOWrapper(handle, encoding="utf-8"), sep="\t")

    expr, summary = hpa_to_expression_long(df, args.source_label, Path(args.url).name)
    expr.to_csv(out, sep="\t", index=False)
    summary.to_csv(out.with_suffix(".summary.tsv"), sep="\t", index=False)
    print(f"Wrote {len(expr):,} HPA expression rows to {out}")


if __name__ == "__main__":
    main()
