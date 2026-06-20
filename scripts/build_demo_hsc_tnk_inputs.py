#!/usr/bin/env python3
"""Build a small HSC/T-NK promoter-design demo dataset.

This is not biological evidence. It is an executable stand-in that validates the
full Stage 2 computational plumbing when curated public hematopoietic activity
tables are not ready yet. Replace it with real promoter_windows.tsv and
activity_long.tsv for results.
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

import pandas as pd

BASES = "ACGT"
CELL_TYPES = ["HSC", "HSPC", "T", "NK", "B", "MYELOID", "ERYTHROID", "MEGAKARYOCYTE"]
ASSAYS = ["accessibility", "initiation", "expression"]
TNK_MOTIFS = ["GGAA", "AGAA", "CTGTGG", "TGAC"]
HSC_MOTIFS = ["GATA", "TGTGGT", "CACCTG"]
OFF_MOTIFS = ["ATGCAA", "TTTGCAT", "CAGCTG"]


def insert_motifs(rng: random.Random, motifs: list[str], length: int = 200) -> str:
    seq = [rng.choice(BASES) for _ in range(length)]
    for motif in motifs:
        pos = rng.randint(10, length - len(motif) - 10)
        seq[pos : pos + len(motif)] = list(motif)
    return "".join(seq)


def activity(kind: str, cell: str, assay: str, rng: random.Random) -> float:
    noise = rng.uniform(0, 0.25)
    if kind == "tnk":
        base = 5.0 if cell in {"T", "NK"} else 0.25
    elif kind == "hsc":
        base = 5.0 if cell in {"HSC", "HSPC"} else 0.4
    elif kind == "offtarget":
        base = 4.0 if cell in {"B", "MYELOID", "ERYTHROID", "MEGAKARYOCYTE"} else 0.4
    else:
        base = 0.2
    if assay == "initiation":
        base *= 0.85
    elif assay == "expression":
        base *= 0.65
    return round(base + noise, 4)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", required=True)
    p.add_argument("--n-per-class", type=int, default=250)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    promoter_rows = []
    activity_rows = []
    classes = [
        ("tnk", TNK_MOTIFS),
        ("hsc", HSC_MOTIFS),
        ("offtarget", OFF_MOTIFS),
        ("inactive", []),
    ]
    idx = 1
    for kind, motifs in classes:
        for _ in range(args.n_per_class):
            pid = f"demo_prom_{idx:06d}"
            chrom = f"chr{rng.randint(1, 22)}"
            start = rng.randint(1_000_000, 200_000_000)
            seq = insert_motifs(rng, motifs)
            promoter_rows.append(
                {
                    "promoter_id": pid,
                    "chr": chrom,
                    "start": start,
                    "end": start + 200,
                    "strand": "+",
                    "tss0": start + 150,
                    "gene_id": f"DEMO{idx:06d}",
                    "gene_name": f"DEMO_{kind}_{idx:06d}",
                    "source": "demo",
                    "sequence": seq,
                    "demo_class": kind,
                }
            )
            for cell in CELL_TYPES:
                for assay in ASSAYS:
                    activity_rows.append(
                        {
                            "promoter_id": pid,
                            "cell_type": cell,
                            "assay": assay,
                            "value": activity(kind, cell, assay, rng),
                        }
                    )
            idx += 1
    pd.DataFrame(promoter_rows).to_csv(out / "promoter_windows.tsv", sep="\t", index=False)
    pd.DataFrame(activity_rows).to_csv(out / "activity_long.tsv", sep="\t", index=False)
    print(f"Wrote demo Stage 2 inputs to {out}")


if __name__ == "__main__":
    main()
