#!/usr/bin/env python3
"""Design an MPRA library with candidate promoters, controls, and barcodes."""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import pandas as pd

from project_utils import read_sequence_file, validate_seq


BASES = "ACGT"


def random_barcode(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def make_barcodes(n: int, length: int, min_distance: int, seed: int) -> list[str]:
    rng = random.Random(seed)
    barcodes: list[str] = []
    attempts = 0
    while len(barcodes) < n:
        attempts += 1
        if attempts > n * 10000:
            raise RuntimeError("Could not make enough barcodes; reduce min distance or lengthen barcode.")
        bc = random_barcode(length, rng)
        if "AAAA" in bc or "TTTT" in bc or "CCCC" in bc or "GGGG" in bc:
            continue
        if all(hamming(bc, other) >= min_distance for other in barcodes):
            barcodes.append(bc)
    return barcodes


def load_controls(paths: list[str]) -> pd.DataFrame:
    rows = []
    for path in paths:
        p = Path(path)
        for i, seq in enumerate(read_sequence_file(p), 1):
            if validate_seq(seq):
                rows.append(
                    {
                        "candidate_id": f"CTRL_{p.stem}_{i:04d}",
                        "sequence": seq,
                        "class": "control",
                        "source_file": p.name,
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranked-candidates", required=True)
    parser.add_argument("--control-files", nargs="*", default=[])
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--barcodes-per-sequence", type=int, default=12)
    parser.add_argument("--barcode-length", type=int, default=16)
    parser.add_argument("--min-barcode-distance", type=int, default=3)
    parser.add_argument("--left-adapter", default="")
    parser.add_argument("--right-adapter", default="")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    cand = pd.read_csv(args.ranked_candidates, sep="\t")
    if "candidate_id" not in cand.columns:
        cand.insert(0, "candidate_id", [f"TNK_CAND_{i+1:05d}" for i in range(len(cand))])
    cand["class"] = cand.get("class", "generated_candidate")
    keep = ["candidate_id", "sequence", "class"]
    for col in ["rank_score", "source_file"]:
        if col in cand.columns:
            keep.append(col)
    library = cand[keep].copy()
    if args.control_files:
        controls = load_controls(args.control_files)
        if len(controls):
            library = pd.concat([library, controls], ignore_index=True, sort=False)
    library = library[library["sequence"].map(validate_seq)].drop_duplicates("candidate_id")

    n_barcodes = len(library) * args.barcodes_per_sequence
    barcodes = make_barcodes(n_barcodes, args.barcode_length, args.min_barcode_distance, args.seed)
    rows = []
    b = 0
    for _, row in library.iterrows():
        for rep in range(args.barcodes_per_sequence):
            barcode = barcodes[b]
            b += 1
            rows.append(
                {
                    **row.to_dict(),
                    "barcode": barcode,
                    "barcode_rep": rep + 1,
                    "oligo": f"{args.left_adapter}{row['sequence']}{args.right_adapter}{barcode}",
                }
            )
    out_prefix = Path(args.output_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_tsv = out_prefix.with_suffix(".library.tsv")
    out_fasta = out_prefix.with_suffix(".oligos.fa")
    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)
    with out_fasta.open("w") as f:
        for r in rows:
            f.write(f">{r['candidate_id']}|bc={r['barcode']}|rep={r['barcode_rep']}\n{r['oligo']}\n")
    print(f"Wrote {len(rows):,} barcode oligos to {out_tsv}")


if __name__ == "__main__":
    main()
