#!/usr/bin/env python3
"""Download public HCA/CELLxGENE marrow data and build expression_long.tsv.

The default source is a public CZ CELLxGENE hematopoietic bone-marrow
collection. The script discovers the H5AD URL through the public collection API,
downloads it, and then reuses the AnnData pseudobulk converter used by Stage 2.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import requests

from project_utils import ensure_dir, write_json


CELLXGENE_COLLECTION_API = "https://api.cellxgene.cziscience.com/curation/v1/collections/{collection_id}"
DEFAULT_COLLECTION_ID = "f6c50495-3361-40ed-a819-fb9644396ed9"


def download_file(url: str, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and out.stat().st_size > 0:
        return out
    with requests.get(url, stream=True, timeout=180) as response:
        response.raise_for_status()
        with out.open("wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                if chunk:
                    handle.write(chunk)
    return out


def fetch_collection(collection_id: str) -> dict:
    url = CELLXGENE_COLLECTION_API.format(collection_id=collection_id)
    response = requests.get(url, timeout=90)
    response.raise_for_status()
    return response.json()


def choose_h5ad_asset(collection: dict, dataset_id: str | None = None) -> tuple[dict, dict]:
    datasets = collection.get("datasets", [])
    if dataset_id:
        datasets = [d for d in datasets if d.get("dataset_id") == dataset_id]
        if not datasets:
            raise SystemExit(f"Dataset id not found in collection: {dataset_id}")
    candidates: list[tuple[int, dict, dict]] = []
    for dataset in datasets:
        title = str(dataset.get("title", "")).lower()
        tissue = " ".join(t.get("label", "") for t in dataset.get("tissue", [])).lower()
        disease = " ".join(d.get("label", "") for d in dataset.get("disease", [])).lower()
        score = 0
        if "bone marrow" in title or "bone marrow" in tissue:
            score += 10
        if "hematopoietic" in title:
            score += 5
        if "normal" in disease:
            score += 2
        score += int(dataset.get("cell_count") or 0) // 100000
        for asset in dataset.get("assets", []):
            if str(asset.get("filetype", "")).upper() == "H5AD" and asset.get("url"):
                candidates.append((score, dataset, asset))
    if not candidates:
        raise SystemExit("No public H5AD asset found in selected CELLxGENE collection.")
    candidates.sort(key=lambda item: (item[0], int(item[2].get("filesize") or 0)), reverse=True)
    _score, dataset, asset = candidates[0]
    return dataset, asset


def run_expression_converter(args: argparse.Namespace, h5ad_path: Path, dataset: dict, collection: dict) -> None:
    converter = Path(__file__).resolve().parent / "expression_from_anndata.py"
    cmd = [
        sys.executable,
        str(converter),
        "--h5ad",
        str(h5ad_path),
        "--output",
        str(args.output),
        "--source",
        args.source_label,
        "--accession",
        dataset.get("dataset_id", ""),
        "--min-cells",
        str(args.min_cells),
        "--seed",
        str(args.seed),
    ]
    if args.cell_type_column:
        cmd.extend(["--cell-type-column", args.cell_type_column])
    if args.gene_id_column:
        cmd.extend(["--gene-id-column", args.gene_id_column])
    if args.gene_name_column:
        cmd.extend(["--gene-name-column", args.gene_name_column])
    if args.layer:
        cmd.extend(["--layer", args.layer])
    if args.use_raw:
        cmd.append("--use-raw")
    if args.max_cells_per_state and args.max_cells_per_state > 0:
        cmd.extend(["--max-cells-per-state", str(args.max_cells_per_state)])
    subprocess.run(cmd, check=True)

    meta = {
        "source": args.source_label,
        "collection_id": collection.get("collection_id"),
        "collection_name": collection.get("name"),
        "collection_url": collection.get("collection_url"),
        "dataset_id": dataset.get("dataset_id"),
        "dataset_title": dataset.get("title"),
        "dataset_cell_count": dataset.get("cell_count"),
        "h5ad_path": str(h5ad_path),
        "expression_long": str(args.output),
        "max_cells_per_state": args.max_cells_per_state,
        "min_cells": args.min_cells,
    }
    write_json(Path(args.output).with_suffix(".metadata.json"), meta)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--collection-id", default=DEFAULT_COLLECTION_ID)
    parser.add_argument("--dataset-id", default=None)
    parser.add_argument("--download-dir", default=Path(__file__).resolve().parents[1] / "downloads" / "hca_cellxgene")
    parser.add_argument("--output", default=Path(__file__).resolve().parents[1] / "data" / "hsc_tnk_real" / "hca_expression_long.tsv")
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--cell-type-column", default=None)
    parser.add_argument("--gene-id-column", default=None)
    parser.add_argument("--gene-name-column", default=None)
    parser.add_argument("--layer", default=None)
    parser.add_argument("--use-raw", action="store_true")
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--max-cells-per-state", type=int, default=5000)
    parser.add_argument("--source-label", default="HCA_CELLXGENE_BONE_MARROW")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    download_dir = ensure_dir(args.download_dir)

    collection = fetch_collection(args.collection_id)
    dataset, asset = choose_h5ad_asset(collection, args.dataset_id)
    h5ad_url = asset["url"]
    h5ad_name = Path(h5ad_url.split("?", 1)[0]).name or f"{dataset['dataset_id']}.h5ad"
    h5ad_path = download_dir / h5ad_name

    manifest = {
        "collection_id": collection.get("collection_id"),
        "collection_name": collection.get("name"),
        "collection_url": collection.get("collection_url"),
        "dataset_id": dataset.get("dataset_id"),
        "dataset_title": dataset.get("title"),
        "dataset_cell_count": dataset.get("cell_count"),
        "asset_filetype": asset.get("filetype"),
        "asset_filesize": asset.get("filesize"),
        "asset_url": h5ad_url,
        "local_h5ad": str(h5ad_path),
    }
    write_json(download_dir / "selected_hca_cellxgene_asset.json", manifest)
    print(json.dumps(manifest, indent=2))

    if args.metadata_only:
        return
    if args.force_download and h5ad_path.exists():
        h5ad_path.unlink()
    download_file(h5ad_url, h5ad_path)
    run_expression_converter(args, h5ad_path, dataset, collection)


if __name__ == "__main__":
    main()
