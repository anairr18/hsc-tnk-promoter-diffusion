#!/usr/bin/env python3
"""Create expression_long.tsv from public HCA/BLUEPRINT-style AnnData files."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_PATTERNS = {
    "HSC": [r"\bhsc\b", r"hematopoietic stem"],
    "HSPC": [r"\bhspc\b", r"progenitor", r"\bmpp\b", r"cd34"],
    "T": [r"\bt cell\b", r"\bt-cell\b", r"\bcd4\b", r"\bcd8\b"],
    "NK": [r"\bnk\b", r"natural killer"],
    "B": [r"\bb cell\b", r"\bb-cell\b", r"naive b", r"memory b"],
    "MYELOID": [r"monocyte", r"myeloid", r"granulocyte", r"neutrophil", r"macrophage"],
    "ERYTHROID": [r"eryth", r"red blood"],
    "MEGAKARYOCYTE": [r"megakary", r"platelet"],
}


def harmonize_label(label: object, patterns: dict[str, list[str]]) -> str | None:
    text = str(label).lower()
    for cell, regexes in patterns.items():
        if any(re.search(pattern, text) for pattern in regexes):
            return cell
    return None


def choose_obs_column(obs: pd.DataFrame, preferred: str | None) -> str:
    if preferred:
        if preferred not in obs.columns:
            raise SystemExit(f"Cell type column not found in h5ad obs: {preferred}")
        return preferred
    candidates = [
        "cell_type",
        "celltype",
        "cell_type__ontology_label",
        "cell_ontology_class",
        "annotation",
        "cluster",
        "leiden",
        "broad_cell_type",
        "lineage",
    ]
    for col in candidates:
        if col in obs.columns:
            return col
    raise SystemExit(f"Could not infer cell type column. Available obs columns include: {list(obs.columns)[:30]}")


def choose_var_column(var: pd.DataFrame, preferred: str | None, fallbacks: list[str]) -> str | None:
    if preferred:
        if preferred not in var.columns:
            raise SystemExit(f"Gene column not found in h5ad var: {preferred}")
        return preferred
    for col in fallbacks:
        if col in var.columns:
            return col
    return None


def matrix_mean(X, indices: np.ndarray) -> np.ndarray:
    sub = X[indices]
    mean = sub.mean(axis=0)
    if hasattr(mean, "A1"):
        return mean.A1
    if hasattr(mean, "A"):
        return np.asarray(mean.A).ravel()
    return np.asarray(mean).ravel()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cell-type-column", default=None)
    parser.add_argument("--gene-id-column", default=None)
    parser.add_argument("--gene-name-column", default=None)
    parser.add_argument("--layer", default=None)
    parser.add_argument("--use-raw", action="store_true")
    parser.add_argument("--cell-map-json", default=None, help="Optional JSON mapping harmonized cell type to regex list.")
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--max-cells-per-state", type=int, default=None)
    parser.add_argument("--source", default="AnnData")
    parser.add_argument("--accession", default="")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    try:
        import anndata as ad
    except ImportError as exc:
        raise SystemExit("Install anndata first: python -m pip install anndata") from exc

    rng = np.random.default_rng(args.seed)
    patterns = DEFAULT_PATTERNS
    if args.cell_map_json:
        patterns = json.loads(Path(args.cell_map_json).read_text())

    adata = ad.read_h5ad(args.h5ad)
    if args.use_raw and adata.raw is not None:
        X = adata.raw.X
        var = adata.raw.var
    elif args.layer:
        X = adata.layers[args.layer]
        var = adata.var
    else:
        X = adata.X
        var = adata.var

    cell_col = choose_obs_column(adata.obs, args.cell_type_column)
    gene_id_col = choose_var_column(var, args.gene_id_column, ["gene_id", "ensembl_id", "feature_id", "id"])
    gene_name_col = choose_var_column(var, args.gene_name_column, ["gene_name", "symbol", "feature_name", "name"])
    gene_ids = var[gene_id_col].astype(str).tolist() if gene_id_col else var.index.astype(str).tolist()
    gene_names = var[gene_name_col].astype(str).tolist() if gene_name_col else var.index.astype(str).tolist()

    labels = adata.obs[cell_col].map(lambda x: harmonize_label(x, patterns))
    rows = []
    summary = []
    for cell_type in sorted(set(labels.dropna())):
        indices = np.where(labels.values == cell_type)[0]
        if len(indices) < args.min_cells:
            continue
        if args.max_cells_per_state and len(indices) > args.max_cells_per_state:
            indices = rng.choice(indices, size=args.max_cells_per_state, replace=False)
        means = matrix_mean(X, indices)
        summary.append({"cell_type": cell_type, "cells": int(len(indices))})
        for gene_id, gene_name, value in zip(gene_ids, gene_names, means):
            rows.append(
                {
                    "cell_type": cell_type,
                    "gene_id": str(gene_id).split(".")[0],
                    "gene_name": gene_name,
                    "value": float(value),
                    "source": args.source,
                    "accession": args.accession,
                    "replicate": f"{cell_type}_pseudobulk",
                }
            )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    pd.DataFrame(summary).to_csv(out.with_suffix(".summary.tsv"), sep="\t", index=False)
    print(f"Wrote {len(rows):,} expression rows to {out}")


if __name__ == "__main__":
    main()
