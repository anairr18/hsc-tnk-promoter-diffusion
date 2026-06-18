#!/usr/bin/env python3
"""Select T/NK-high, HSC/non-target-low promoter training candidates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from project_utils import validate_seq


DEFAULT_TARGETS = ["T", "CD4_T", "CD8_T", "NK"]
DEFAULT_OFFTARGETS = [
    "HSC",
    "LT_HSC",
    "HSPC",
    "MPP",
    "B",
    "MYELOID",
    "MONOCYTE",
    "GRANULOCYTE",
    "ERYTHROID",
    "MEGAKARYOCYTE",
]


def split_csv(value: str | None, default: list[str]) -> list[str]:
    if not value:
        return default
    return [x.strip() for x in value.split(",") if x.strip()]


def feature_cols_for_cells(columns: list[str], cells: list[str]) -> list[str]:
    out = []
    lower_cells = [c.lower() for c in cells]
    for col in columns:
        if "__" not in col:
            continue
        assay, cell = col.split("__", 1)
        cell_lower = cell.lower()
        if any(c == cell_lower or c in cell_lower for c in lower_cells):
            out.append(col)
    return out


def robust_z(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df[cols].astype(float).copy()
    for c in cols:
        med = out[c].median()
        mad = (out[c] - med).abs().median()
        denom = mad * 1.4826 if mad > 0 else out[c].std()
        out[c] = 0.0 if not denom or np.isnan(denom) else (out[c] - med) / denom
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--activity-matrix", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--profile-json", default=None)
    parser.add_argument("--target-cells", default=None)
    parser.add_argument("--offtarget-cells", default=None)
    parser.add_argument("--top-n", type=int, default=5000)
    parser.add_argument("--min-specificity", type=float, default=0.0)
    args = parser.parse_args()

    df = pd.read_csv(args.activity_matrix, sep="\t")
    if "sequence" not in df.columns:
        raise SystemExit("activity matrix must include a sequence column")
    df = df[df["sequence"].map(validate_seq)].copy()

    target_cells = split_csv(args.target_cells, DEFAULT_TARGETS)
    offtarget_cells = split_csv(args.offtarget_cells, DEFAULT_OFFTARGETS)
    all_feature_cols = [c for c in df.columns if "__" in c]
    target_cols = feature_cols_for_cells(all_feature_cols, target_cells)
    off_cols = feature_cols_for_cells(all_feature_cols, offtarget_cells)
    if not target_cols:
        raise SystemExit(f"No target feature columns found for cells: {target_cells}")
    if not off_cols:
        raise SystemExit(f"No off-target feature columns found for cells: {offtarget_cells}")

    z = robust_z(df, all_feature_cols)
    df["tnk_target_score"] = z[target_cols].mean(axis=1)
    df["offtarget_max_score"] = z[off_cols].max(axis=1)
    df["tnk_specificity_score"] = df["tnk_target_score"] - df["offtarget_max_score"]
    df["TAG"] = "TNK_HIGH"

    selected = df[df["tnk_specificity_score"] >= args.min_specificity].copy()
    selected = selected.sort_values("tnk_specificity_score", ascending=False).head(args.top_n)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(out, sep="\t", index=False)

    profile_path = Path(args.profile_json) if args.profile_json else out.with_suffix(".profiles.json")
    profiles = {
        "target_cells": target_cells,
        "offtarget_cells": offtarget_cells,
        "profile_columns": all_feature_cols,
        "target_columns": target_cols,
        "offtarget_columns": off_cols,
        "selected_n": int(len(selected)),
    }
    profile_path.write_text(json.dumps(profiles, indent=2) + "\n")
    print(f"Wrote {len(selected):,} selected promoters to {out}")


if __name__ == "__main__":
    main()
