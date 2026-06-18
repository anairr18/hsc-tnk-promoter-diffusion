#!/usr/bin/env python3
"""Extract strand-aware 200bp promoter-proximal windows."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd

from project_utils import revcomp, validate_seq


def parse_attrs(attr_text: str) -> dict[str, str]:
    out = {}
    for item in attr_text.split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        out[key] = value.strip().strip('"')
    return out


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def load_gencode_tss(path: Path) -> pd.DataFrame:
    rows = []
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parse_attrs(parts[8])
            gene_type = attrs.get("gene_type") or attrs.get("gene_biotype")
            if gene_type and gene_type != "protein_coding":
                continue
            chrom = parts[0] if parts[0].startswith("chr") else f"chr{parts[0]}"
            start_1 = int(parts[3])
            end_1 = int(parts[4])
            strand = parts[6]
            tss0 = start_1 - 1 if strand == "+" else end_1 - 1
            rows.append(
                {
                    "chr": chrom,
                    "tss0": tss0,
                    "strand": strand,
                    "gene_id": attrs.get("gene_id", "").split(".")[0],
                    "gene_name": attrs.get("gene_name", ""),
                    "source": "GENCODE",
                }
            )
    return pd.DataFrame(rows).drop_duplicates(["chr", "tss0", "strand", "gene_id"])


def load_tss_bed(path: Path, source: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 6:
        raise ValueError("TSS BED must have at least 6 columns: chr start end name score strand")
    rows = []
    for _, row in df.iterrows():
        chrom, start, end, name, _, strand = row.iloc[:6]
        tss0 = int(start) if strand == "+" else int(end) - 1
        rows.append(
            {
                "chr": chrom if str(chrom).startswith("chr") else f"chr{chrom}",
                "tss0": tss0,
                "strand": strand,
                "gene_id": str(name),
                "gene_name": str(name),
                "source": source,
            }
        )
    return pd.DataFrame(rows)


def extract_sequence(fasta, chrom: str, start: int, end: int, strand: str) -> str | None:
    try:
        seq = str(fasta[chrom][start:end]).upper()
    except Exception:
        return None
    if strand == "-":
        seq = revcomp(seq)
    return seq if validate_seq(seq, end - start) else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, help="hg38 FASTA indexed by pyfaidx.")
    parser.add_argument("--gencode-gtf", default=None)
    parser.add_argument("--fantom-tss-bed", default=None)
    parser.add_argument("--output", required=True)
    parser.add_argument("--upstream", type=int, default=150)
    parser.add_argument("--downstream", type=int, default=50)
    args = parser.parse_args()

    if not args.gencode_gtf and not args.fantom_tss_bed:
        raise SystemExit("Provide --fantom-tss-bed, --gencode-gtf, or both.")
    try:
        from pyfaidx import Fasta
    except ImportError as exc:
        raise SystemExit("Install pyfaidx first: python -m pip install pyfaidx") from exc

    frames = []
    if args.fantom_tss_bed:
        frames.append(load_tss_bed(Path(args.fantom_tss_bed), "FANTOM5"))
    if args.gencode_gtf:
        frames.append(load_gencode_tss(Path(args.gencode_gtf)))
    tss = pd.concat(frames, ignore_index=True).drop_duplicates(["chr", "tss0", "strand", "gene_id"])

    fasta = Fasta(args.fasta)
    rows = []
    for i, row in tss.iterrows():
        if row["strand"] == "+":
            start = max(0, int(row["tss0"]) - args.upstream)
            end = int(row["tss0"]) + args.downstream
        else:
            start = max(0, int(row["tss0"]) - args.downstream)
            end = int(row["tss0"]) + args.upstream
        seq = extract_sequence(fasta, row["chr"], start, end, row["strand"])
        if not seq:
            continue
        rows.append(
            {
                "promoter_id": f"prom_{len(rows)+1:08d}",
                "chr": row["chr"],
                "start": start,
                "end": end,
                "strand": row["strand"],
                "tss0": int(row["tss0"]),
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
                "source": row["source"],
                "sequence": seq,
            }
        )
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"Wrote {len(rows):,} promoter windows to {out}")


if __name__ == "__main__":
    main()
