#!/usr/bin/env python3
"""Merge promoter windows with accessibility/initiation/expression activity."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REQUIRED_LONG_COLUMNS = {"promoter_id", "cell_type", "assay", "value"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--promoters", required=True, help="TSV from build_promoter_windows.py")
    parser.add_argument("--activity-long", required=True, help="Long TSV: promoter_id, cell_type, assay, value")
    parser.add_argument("--output", required=True)
    parser.add_argument("--fill-value", type=float, default=0.0)
    args = parser.parse_args()

    promoters = pd.read_csv(args.promoters, sep="\t")
    activity = pd.read_csv(args.activity_long, sep="\t")
    missing = REQUIRED_LONG_COLUMNS - set(activity.columns)
    if missing:
        raise SystemExit(f"--activity-long missing columns: {sorted(missing)}")

    activity = activity.copy()
    activity["feature"] = activity["assay"].astype(str) + "__" + activity["cell_type"].astype(str)
    wide = activity.pivot_table(
        index="promoter_id",
        columns="feature",
        values="value",
        aggfunc="median",
    ).reset_index()
    wide.columns.name = None
    merged = promoters.merge(wide, on="promoter_id", how="left")
    feature_cols = [c for c in merged.columns if "__" in c]
    merged[feature_cols] = merged[feature_cols].fillna(args.fill_value)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, sep="\t", index=False)
    print(f"Wrote {len(merged):,} promoters x {len(feature_cols):,} activity features to {out}")


if __name__ == "__main__":
    main()
