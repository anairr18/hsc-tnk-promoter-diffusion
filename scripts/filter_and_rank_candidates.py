#!/usr/bin/env python3
"""Filter and rank generated promoter candidates for MPRA nomination."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from project_utils import cpg_count, gc_content, read_sequence_file, validate_seq


UNWANTED_MOTIFS = ["AAAAAA", "TTTTTT", "GCGCGCGC"]


def read_candidates(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".tsv", ".csv"}:
        sep = "\t" if path.suffix.lower() == ".tsv" else ","
        df = pd.read_csv(path, sep=sep)
        if "sequence" in df.columns:
            return df
    seqs = read_sequence_file(path)
    return pd.DataFrame({"sequence": seqs})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True, nargs="+")
    parser.add_argument("--predictions", default=None, help="Optional TSV with sequence and prediction columns.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-gc", type=float, default=0.30)
    parser.add_argument("--max-gc", type=float, default=0.75)
    parser.add_argument("--max-cpg", type=int, default=35)
    parser.add_argument("--top-n", type=int, default=600)
    args = parser.parse_args()

    frames = []
    for item in args.candidates:
        df = read_candidates(Path(item))
        df["source_file"] = Path(item).name
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df["sequence"] = df["sequence"].astype(str).str.upper().str.strip()
    df = df[df["sequence"].map(validate_seq)].drop_duplicates("sequence").copy()
    df["gc"] = df["sequence"].map(gc_content)
    df["cpg"] = df["sequence"].map(cpg_count)
    df["has_unwanted_motif"] = df["sequence"].map(lambda s: any(m in s for m in UNWANTED_MOTIFS))
    df = df[
        (df["gc"] >= args.min_gc)
        & (df["gc"] <= args.max_gc)
        & (df["cpg"] <= args.max_cpg)
        & (~df["has_unwanted_motif"])
    ].copy()

    if args.predictions:
        pred = pd.read_csv(args.predictions, sep="\t")
        df = df.merge(pred, on="sequence", how="left", suffixes=("", "_pred"))

    target_cols = [c for c in df.columns if "tnk" in c.lower() and "score" in c.lower()]
    off_cols = [c for c in df.columns if any(x in c.lower() for x in ["hsc", "offtarget", "myeloid", "bcell", "erythroid"])]
    if target_cols:
        df["rank_score"] = df[target_cols].mean(axis=1)
        if off_cols:
            df["rank_score"] = df["rank_score"] - df[off_cols].max(axis=1)
    else:
        # Heuristic fallback: prefer moderate GC/CpG and motif-safe novelty.
        df["rank_score"] = 1.0 - (df["gc"] - 0.52).abs() - (df["cpg"] / 100.0)

    df = df.sort_values("rank_score", ascending=False).head(args.top_n)
    df.insert(0, "candidate_id", [f"TNK_CAND_{i+1:05d}" for i in range(len(df))])
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"Wrote {len(df):,} ranked candidates to {out}")


if __name__ == "__main__":
    main()
