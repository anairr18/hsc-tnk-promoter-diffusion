#!/usr/bin/env python3
"""Create positive, leakiness, inactive, and scrambled MPRA controls."""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import pandas as pd


def shuffle_seq(seq: str, rng: random.Random) -> str:
    letters = list(seq)
    rng.shuffle(letters)
    return "".join(letters)


def feature_cols(df: pd.DataFrame, assays: list[str], cells: list[str]) -> list[str]:
    cols = []
    for c in df.columns:
        low = c.lower()
        if "__" not in c:
            continue
        if any(a.lower() in low for a in assays) and any(cell.lower() in low for cell in cells):
            cols.append(c)
    return cols


def write_txt(path: Path, seqs: list[str]) -> None:
    path.write_text("\n".join(seqs) + ("\n" if seqs else ""))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--activity-matrix", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--n", type=int, default=96)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    rng = random.Random(args.seed)
    df = pd.read_csv(args.activity_matrix, sep="\t")
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    tnk_cols = feature_cols(df, ["accessibility", "initiation", "expression"], ["T", "NK"])
    hsc_cols = feature_cols(df, ["accessibility", "initiation", "expression"], ["HSC", "HSPC"])
    off_cols = feature_cols(df, ["accessibility", "initiation", "expression"], ["B", "MYELOID", "ERYTHROID", "MEGAKARYOCYTE"])
    all_cols = [c for c in df.columns if "__" in c]

    df = df.copy()
    df["tnk_control_score"] = df[tnk_cols].mean(axis=1) if tnk_cols else 0.0
    df["hsc_control_score"] = df[hsc_cols].mean(axis=1) if hsc_cols else 0.0
    df["offtarget_control_score"] = df[off_cols].mean(axis=1) if off_cols else 0.0
    df["inactive_score"] = df[all_cols].mean(axis=1) if all_cols else 0.0

    positives = df.sort_values("tnk_control_score", ascending=False)["sequence"].head(args.n).tolist()
    hsc_leaky = df.sort_values("hsc_control_score", ascending=False)["sequence"].head(args.n).tolist()
    inactive = df.sort_values("inactive_score", ascending=True)["sequence"].head(args.n).tolist()
    scrambled = [shuffle_seq(s, rng) for s in positives[: args.n]]

    write_txt(out / "tnk_positive_controls.txt", positives)
    write_txt(out / "hsc_leakiness_controls.txt", hsc_leaky)
    write_txt(out / "inactive_controls.txt", inactive)
    write_txt(out / "scrambled_controls.txt", scrambled)
    print(f"Wrote MPRA control files to {out}")


if __name__ == "__main__":
    main()
