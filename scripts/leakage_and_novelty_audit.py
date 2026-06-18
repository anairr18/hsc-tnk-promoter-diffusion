#!/usr/bin/env python3
"""Audit exact and approximate leakage between train/test/generated sequences."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from project_utils import kmer_counts, normalize_seq, read_sequence_file, validate_seq


def kmerset(seq: str, k: int) -> set[str]:
    seq = normalize_seq(seq)
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def max_jaccard(query: str, refs: list[set[str]], k: int) -> float:
    q = kmerset(query, k)
    if not q or not refs:
        return float("nan")
    best = 0.0
    for r in refs:
        union = len(q | r)
        if union:
            best = max(best, len(q & r) / union)
    return best


def load_split(path: Path) -> list[str]:
    df = pd.read_csv(path, sep="\t")
    if "sequence" in df.columns:
        return df["sequence"].astype(str).str.upper().tolist()
    return read_sequence_file(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True)
    parser.add_argument("--val", default=None)
    parser.add_argument("--test", required=True)
    parser.add_argument("--generated", nargs="*", default=[])
    parser.add_argument("--output", required=True)
    parser.add_argument("--k", type=int, default=10)
    args = parser.parse_args()

    train = [normalize_seq(s) for s in load_split(Path(args.train)) if validate_seq(s)]
    val = [normalize_seq(s) for s in load_split(Path(args.val)) if validate_seq(s)] if args.val else []
    test = [normalize_seq(s) for s in load_split(Path(args.test)) if validate_seq(s)]
    train_set = set(train)
    val_set = set(val)
    train_kmers = [kmerset(s, args.k) for s in train]

    rows = []
    for name, seqs in [("test", test), ("val", val)]:
        if not seqs:
            continue
        exact_train = sum(1 for s in seqs if s in train_set)
        exact_other = sum(1 for s in seqs if s in val_set) if name == "test" else 0
        approx = [max_jaccard(s, train_kmers, args.k) for s in seqs[: min(2000, len(seqs))]]
        rows.append(
            {
                "group": name,
                "n": len(seqs),
                "exact_train_overlap": exact_train,
                "exact_train_overlap_rate": exact_train / len(seqs),
                "exact_val_overlap": exact_other,
                "max_train_jaccard_kmer_sample_mean": sum(approx) / len(approx) if approx else float("nan"),
                "max_train_jaccard_kmer_sample_max": max(approx) if approx else float("nan"),
            }
        )
    for gen_path in args.generated:
        seqs = [normalize_seq(s) for s in read_sequence_file(gen_path) if validate_seq(s)]
        if not seqs:
            continue
        exact_train = sum(1 for s in seqs if s in train_set)
        approx = [max_jaccard(s, train_kmers, args.k) for s in seqs[: min(2000, len(seqs))]]
        rows.append(
            {
                "group": Path(gen_path).stem,
                "n": len(seqs),
                "exact_train_overlap": exact_train,
                "exact_train_overlap_rate": exact_train / len(seqs),
                "exact_val_overlap": 0,
                "max_train_jaccard_kmer_sample_mean": sum(approx) / len(approx) if approx else float("nan"),
                "max_train_jaccard_kmer_sample_max": max(approx) if approx else float("nan"),
            }
        )
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"Wrote leakage audit to {out}")


if __name__ == "__main__":
    main()
