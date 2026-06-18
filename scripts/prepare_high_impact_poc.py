#!/usr/bin/env python3
"""Prepare publication-grade POC split bundles for DNA-Diffusion.

This creates multiple benchmark splits from the same 200bp sequence table:
  - random stratified split for continuity with the current POC
  - chromosome-heldout split to test genomic generalization
  - gene-heldout split to reduce promoter/enhancer-neighborhood leakage

Each split is exported as DNA-Diffusion-compatible TSVs and a saved-data pickle
so the upstream dataloader uses the intended split instead of re-splitting by
chr1/chr2 internally.
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import pandas as pd
from sklearn.model_selection import GroupShuffleSplit, train_test_split

from project_utils import ensure_dir, validate_seq, write_json

REQUIRED = {"chr", "sequence", "TAG"}


def normalize_table(df: pd.DataFrame) -> pd.DataFrame:
    missing = REQUIRED - set(df.columns)
    if missing:
        raise SystemExit(f"Input table missing required columns: {sorted(missing)}")
    df = df.copy()
    df["sequence"] = df["sequence"].astype(str).str.upper().str.strip()
    df = df[df["sequence"].map(validate_seq)].drop_duplicates(["sequence", "TAG"]).reset_index(drop=True)
    if "gene" not in df.columns:
        if "gene_name" in df.columns:
            df["gene"] = df["gene_name"].fillna("")
        elif "gene_id" in df.columns:
            df["gene"] = df["gene_id"].fillna("")
        else:
            df["gene"] = df["chr"].astype(str) + ":" + (df.index // 5).astype(str)
    return df


def random_split(df: pd.DataFrame, seed: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    train_val, test = train_test_split(df, test_size=0.1, random_state=seed, stratify=df["TAG"])
    train, val = train_test_split(train_val, test_size=0.111111, random_state=seed, stratify=train_val["TAG"])
    return train, val, test


def chromosome_split(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    test_chr = {"chr1", "chr8", "chr21"}
    val_chr = {"chr2", "chr9", "chr22"}
    test = df[df["chr"].isin(test_chr)]
    val = df[df["chr"].isin(val_chr)]
    train = df[~df["chr"].isin(test_chr | val_chr)]
    if min(len(train), len(val), len(test)) == 0:
        raise ValueError("Chromosome split produced an empty split; check chr names.")
    return train, val, test


def group_split(df: pd.DataFrame, seed: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    groups = df["gene"].astype(str)
    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed)
    train_idx, rest_idx = next(gss.split(df, df["TAG"], groups=groups))
    train = df.iloc[train_idx]
    rest = df.iloc[rest_idx]
    rest_groups = rest["gene"].astype(str)
    gss2 = GroupShuffleSplit(n_splits=1, test_size=0.5, random_state=seed)
    val_idx, test_idx = next(gss2.split(rest, rest["TAG"], groups=rest_groups))
    return train, rest.iloc[val_idx], rest.iloc[test_idx]


def export_split(name: str, split_dir: Path, train: pd.DataFrame, val: pd.DataFrame, test: pd.DataFrame) -> dict:
    out = ensure_dir(split_dir / name)
    for split_name, split_df in [("train", train), ("val", val), ("test", test)]:
        split_df[["chr", "sequence", "TAG"]].to_csv(out / f"{split_name}.txt", sep="\t", index=False)
        split_df.to_csv(out / f"{split_name}_with_metadata.tsv", sep="\t", index=False)
    with (out / "encode_data.pkl").open("wb") as f:
        pickle.dump(
            {
                "train_df": train[["chr", "sequence", "TAG"]].reset_index(drop=True),
                "validation_df": val[["chr", "sequence", "TAG"]].reset_index(drop=True),
                "test_df": test[["chr", "sequence", "TAG"]].reset_index(drop=True),
            },
            f,
        )
    summary = {
        "split": name,
        "train_n": int(len(train)),
        "val_n": int(len(val)),
        "test_n": int(len(test)),
        "train_labels": train["TAG"].value_counts().to_dict(),
        "val_labels": val["TAG"].value_counts().to_dict(),
        "test_labels": test["TAG"].value_counts().to_dict(),
    }
    write_json(out / "split_summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Full sequence table with chr, sequence, TAG.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    input_path = Path(args.input)
    sep = "\t" if input_path.suffix.lower() in {".tsv", ".txt"} else ","
    df = normalize_table(pd.read_csv(input_path, sep=sep))
    out = ensure_dir(args.output_dir)

    summaries = []
    train, val, test = random_split(df, args.seed)
    summaries.append(export_split("random_stratified", out, train, val, test))
    train, val, test = chromosome_split(df)
    summaries.append(export_split("chromosome_holdout", out, train, val, test))
    train, val, test = group_split(df, args.seed)
    summaries.append(export_split("gene_holdout", out, train, val, test))

    write_json(out / "all_split_summaries.json", summaries)
    print(f"Wrote {len(summaries)} split bundles to {out}")


if __name__ == "__main__":
    main()
