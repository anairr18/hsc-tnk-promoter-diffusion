#!/usr/bin/env python3
"""Create the next design round from MPRA activity measurements."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mpra-hits", required=True, help="sequence_level_hits.tsv from analyze_mpra_barcodes.py")
    parser.add_argument("--candidate-features", default=None)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--top-positive", type=int, default=200)
    parser.add_argument("--top-negative", type=int, default=200)
    args = parser.parse_args()

    hits = pd.read_csv(args.mpra_hits, sep="\t")
    required = {"candidate_id", "sequence", "tnk_specificity", "hsc_leakiness_delta"}
    missing = required - set(hits.columns)
    if missing:
        raise SystemExit(f"MPRA hits table missing columns: {sorted(missing)}")
    hits["active_learning_score"] = hits["tnk_specificity"].fillna(-999) + hits["hsc_leakiness_delta"].fillna(-999)
    positives = hits.sort_values("active_learning_score", ascending=False).head(args.top_positive).copy()
    negatives = hits.sort_values("active_learning_score", ascending=True).head(args.top_negative).copy()
    positives["round_label"] = "positive_tnk_specific"
    negatives["round_label"] = "negative_or_leaky"
    training = pd.concat([positives, negatives], ignore_index=True)

    if args.candidate_features:
        feats = pd.read_csv(args.candidate_features, sep="\t")
        shared = [c for c in ["candidate_id", "sequence"] if c in feats.columns and c in training.columns]
        training = training.merge(feats, on=shared, how="left", suffixes=("", "_feature"))

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    positives.to_csv(out / "round2_positive_seeds.tsv", sep="\t", index=False)
    negatives.to_csv(out / "round2_negative_seeds.tsv", sep="\t", index=False)
    training.to_csv(out / "round2_active_learning_training_table.tsv", sep="\t", index=False)
    print(f"Wrote active-learning round files to {out}")


if __name__ == "__main__":
    main()
