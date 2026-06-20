#!/usr/bin/env python3
"""Annotate promoter windows with overlapping ENCODE SCREEN cCRE classes."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


def load_ccre_bed(path: Path, class_column: int, accession_column: int) -> dict[str, list[tuple[int, int, str, str]]]:
    intervals: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            cls = parts[class_column - 1] if len(parts) >= class_column else parts[3] if len(parts) >= 4 else "cCRE"
            accession = parts[accession_column - 1] if len(parts) >= accession_column else parts[3] if len(parts) >= 4 else ""
            intervals[chrom].append((start, end, accession, cls))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start < b_end and b_start < a_end


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--promoters", required=True)
    parser.add_argument("--ccre-bed", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--class-column", type=int, default=10, help="1-based BED column containing cCRE class; fallback to name.")
    parser.add_argument("--accession-column", type=int, default=4, help="1-based BED column containing cCRE accession/name.")
    args = parser.parse_args()

    promoters = pd.read_csv(args.promoters, sep="\t")
    required = {"chr", "start", "end"}
    missing = required - set(promoters.columns)
    if missing:
        raise SystemExit(f"promoters missing columns: {sorted(missing)}")
    intervals = load_ccre_bed(Path(args.ccre_bed), args.class_column, args.accession_column)
    rows = []
    for _, prom in promoters.iterrows():
        chrom = str(prom["chr"])
        p_start = int(prom["start"])
        p_end = int(prom["end"])
        hits = []
        for c_start, c_end, accession, cls in intervals.get(chrom, []):
            if c_start >= p_end:
                break
            if overlaps(p_start, p_end, c_start, c_end):
                hits.append((accession, cls))
        row = prom.to_dict()
        row["ccre_accessions"] = ";".join(sorted({h[0] for h in hits if h[0]}))
        row["ccre_classes"] = ";".join(sorted({h[1] for h in hits if h[1]}))
        row["ccre_overlap_count"] = len(hits)
        rows.append(row)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"Wrote cCRE annotations for {len(rows):,} promoters to {out}")


if __name__ == "__main__":
    main()
