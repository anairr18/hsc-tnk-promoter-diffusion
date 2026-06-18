#!/usr/bin/env python3
"""Analyze MPRA RNA/DNA barcode counts for HSC-to-T/NK validation."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COUNTS = {"barcode", "sample_id", "donor", "timepoint", "cell_type", "molecule", "count"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True, help="Long TSV of barcode counts.")
    parser.add_argument("--library", required=True, help="MPRA library TSV from design_mpra_library.py.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--target-cells", default="T,NK")
    parser.add_argument("--hsc-labels", default="CD34,HSC,HSPC")
    parser.add_argument("--pseudocount", type=float, default=1.0)
    args = parser.parse_args()

    counts = pd.read_csv(args.counts, sep="\t")
    missing = REQUIRED_COUNTS - set(counts.columns)
    if missing:
        raise SystemExit(f"--counts missing columns: {sorted(missing)}")
    library = pd.read_csv(args.library, sep="\t")
    counts = counts.merge(library[["barcode", "candidate_id", "sequence", "class"]], on="barcode", how="inner")

    pivot = counts.pivot_table(
        index=["candidate_id", "sequence", "class", "barcode", "sample_id", "donor", "timepoint", "cell_type"],
        columns="molecule",
        values="count",
        aggfunc="sum",
        fill_value=0,
    ).reset_index()
    pivot.columns.name = None
    if "RNA" not in pivot.columns or "DNA" not in pivot.columns:
        raise SystemExit("counts must include molecule values RNA and DNA")
    pivot["log2_rna_dna"] = np.log2((pivot["RNA"] + args.pseudocount) / (pivot["DNA"] + args.pseudocount))

    barcode_summary = pivot.groupby(
        ["candidate_id", "sequence", "class", "donor", "timepoint", "cell_type"], as_index=False
    )["log2_rna_dna"].agg(["mean", "std", "count"]).reset_index()
    barcode_summary = barcode_summary.rename(columns={"mean": "activity_log2_rna_dna", "std": "activity_sd", "count": "barcode_n"})

    targets = [x.strip().lower() for x in args.target_cells.split(",") if x.strip()]
    hsc_labels = [x.strip().lower() for x in args.hsc_labels.split(",") if x.strip()]
    barcode_summary["is_target"] = barcode_summary["cell_type"].str.lower().map(
        lambda x: any(t in x for t in targets)
    )
    barcode_summary["is_hsc"] = barcode_summary["cell_type"].str.lower().map(
        lambda x: any(h in x for h in hsc_labels)
    )

    seq_rows = []
    for cid, sub in barcode_summary.groupby("candidate_id"):
        target = sub.loc[sub["is_target"], "activity_log2_rna_dna"]
        hsc = sub.loc[sub["is_hsc"], "activity_log2_rna_dna"]
        off = sub.loc[~sub["is_target"], "activity_log2_rna_dna"]
        seq_rows.append(
            {
                "candidate_id": cid,
                "sequence": sub["sequence"].iloc[0],
                "class": sub["class"].iloc[0],
                "target_activity": target.mean() if len(target) else np.nan,
                "hsc_activity": hsc.mean() if len(hsc) else np.nan,
                "offtarget_activity": off.mean() if len(off) else np.nan,
                "tnk_specificity": (target.mean() - off.mean()) if len(target) and len(off) else np.nan,
                "hsc_leakiness_delta": (target.mean() - hsc.mean()) if len(target) and len(hsc) else np.nan,
                "measurements_n": len(sub),
            }
        )
    sequence_summary = pd.DataFrame(seq_rows).sort_values(
        ["tnk_specificity", "hsc_leakiness_delta"], ascending=False, na_position="last"
    )

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    pivot.to_csv(out / "barcode_level_activity.tsv", sep="\t", index=False)
    barcode_summary.to_csv(out / "sample_level_activity.tsv", sep="\t", index=False)
    sequence_summary.to_csv(out / "sequence_level_hits.tsv", sep="\t", index=False)
    print(f"Wrote MPRA analysis tables to {out}")


if __name__ == "__main__":
    main()
