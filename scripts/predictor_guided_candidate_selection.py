#!/usr/bin/env python3
"""Rank candidates with an explicit high-target/low-offtarget reward."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from project_utils import cpg_count, gc_content, read_sequence_file, validate_seq


DEFAULT_WEIGHTS = {
    "target_regex": ["tnk", "t_cell", "nk", "target_activity"],
    "offtarget_regex": ["hsc", "hspc", "b_cell", "myeloid", "erythroid", "offtarget"],
    "target_weight": 1.0,
    "offtarget_weight": 1.0,
    "gc_penalty_weight": 0.5,
    "cpg_penalty_weight": 0.02,
    "gc_center": 0.52,
}


def load_candidates(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".tsv", ".csv"}:
        df = pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else ",")
        if "sequence" in df.columns:
            return df
    return pd.DataFrame({"sequence": read_sequence_file(path)})


def matching_columns(df: pd.DataFrame, patterns: list[str]) -> list[str]:
    out = []
    for c in df.columns:
        low = c.lower()
        if any(p.lower() in low for p in patterns):
            if pd.api.types.is_numeric_dtype(df[c]):
                out.append(c)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--scores", required=True, help="Merged predictor scores TSV.")
    parser.add_argument("--weights-json", default=None)
    parser.add_argument("--output", required=True)
    parser.add_argument("--top-n", type=int, default=600)
    args = parser.parse_args()

    weights = DEFAULT_WEIGHTS.copy()
    if args.weights_json:
        weights.update(json.loads(Path(args.weights_json).read_text()))

    cand = load_candidates(Path(args.candidates))
    cand["sequence"] = cand["sequence"].astype(str).str.upper().str.strip()
    cand = cand[cand["sequence"].map(validate_seq)].drop_duplicates("sequence")
    scores = pd.read_csv(args.scores, sep="\t")
    if "sequence" not in scores.columns:
        raise SystemExit("--scores must include a sequence column")
    df = cand.merge(scores, on="sequence", how="left")
    df["gc"] = df["sequence"].map(gc_content)
    df["cpg"] = df["sequence"].map(cpg_count)

    target_cols = matching_columns(df, weights["target_regex"])
    off_cols = matching_columns(df, weights["offtarget_regex"])
    if not target_cols:
        raise SystemExit(f"No numeric target columns matched: {weights['target_regex']}")
    df["target_score"] = df[target_cols].mean(axis=1)
    df["offtarget_score"] = df[off_cols].max(axis=1) if off_cols else 0.0
    df["gc_penalty"] = (df["gc"] - float(weights["gc_center"])).abs()
    df["cpg_penalty"] = df["cpg"]
    df["guided_reward"] = (
        float(weights["target_weight"]) * df["target_score"]
        - float(weights["offtarget_weight"]) * df["offtarget_score"]
        - float(weights["gc_penalty_weight"]) * df["gc_penalty"]
        - float(weights["cpg_penalty_weight"]) * df["cpg_penalty"]
    )
    df = df.sort_values("guided_reward", ascending=False).head(args.top_n)
    df.insert(0, "guided_rank", range(1, len(df) + 1))
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    (out.with_suffix(".weights.json")).write_text(json.dumps(weights, indent=2) + "\n")
    print(f"Wrote {len(df):,} guided candidates to {out}")


if __name__ == "__main__":
    main()
