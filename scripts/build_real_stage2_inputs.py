#!/usr/bin/env python3
"""Build real HSC/T-NK Stage 2 promoter windows and activity tables.

This script turns public source files into the two required real-data inputs:

- promoter_windows.tsv
- activity_long.tsv

It intentionally accepts already-downloaded public data files rather than
hard-coding every portal-specific download format. References can be downloaded
automatically; omics tracks are provided through simple manifests.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import shutil
from collections import defaultdict
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import requests

from project_utils import ensure_dir, revcomp, validate_seq, write_json


HG38_FASTA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
GENCODE_V38_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz"
VALID_CELL_TYPES = {"HSC", "HSPC", "T", "NK", "B", "MYELOID", "ERYTHROID", "MEGAKARYOCYTE"}
VALID_ASSAYS = {"accessibility", "initiation", "expression"}
REQUIRED_ACTIVITY_COLUMNS = {"promoter_id", "cell_type", "assay", "value"}


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def parse_attrs(attr_text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in attr_text.split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, value = item.split(" ", 1)
        out[key] = value.strip().strip('"')
    return out


def download(url: str, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists() and out.stat().st_size > 0:
        return out
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with out.open("wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                if chunk:
                    handle.write(chunk)
    return out


def maybe_gunzip(path: Path) -> Path:
    if path.suffix != ".gz":
        return path
    out = path.with_suffix("")
    if out.exists() and out.stat().st_size > 0:
        return out
    with gzip.open(path, "rb") as src, out.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return out


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_row(path: Path, role: str, source: str = "", accession: str = "", skip_checksum: bool = False) -> dict[str, str | int]:
    return {
        "role": role,
        "source": source,
        "accession": accession,
        "path": str(path),
        "file_name": path.name,
        "bytes": path.stat().st_size if path.exists() else 0,
        "sha256": "" if skip_checksum or not path.exists() else sha256_file(path),
    }


def stable_promoter_id(chrom: str, start: int, end: int, strand: str, gene_key: str) -> str:
    key = f"{chrom}:{start}-{end}:{strand}:{gene_key}"
    return "prom_" + hashlib.sha1(key.encode("utf-8")).hexdigest()[:14]


def normalize_activity_value(value: object) -> float:
    try:
        val = float(value)
    except Exception:
        return math.nan
    if math.isnan(val):
        return math.nan
    return math.log1p(max(0.0, val))


def harmonize_cell(value: str) -> str:
    cell = str(value).strip().upper().replace("-", "_").replace(" ", "_")
    aliases = {
        "CD34": "HSC",
        "CD34_POSITIVE": "HSC",
        "HSC_HSPC": "HSPC",
        "CD4_T": "T",
        "CD8_T": "T",
        "TCELL": "T",
        "T_CELL": "T",
        "NK_CELL": "NK",
        "BCELL": "B",
        "B_CELL": "B",
        "MONOCYTE": "MYELOID",
        "GRANULOCYTE": "MYELOID",
        "MEGAKARYOCYTE_ERYTHROID": "ERYTHROID",
    }
    return aliases.get(cell, cell)


def harmonize_assay(value: str) -> str:
    assay = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "atac": "accessibility",
        "atac_seq": "accessibility",
        "dnase": "accessibility",
        "dnase_seq": "accessibility",
        "cage": "initiation",
        "rampage": "initiation",
        "tss": "initiation",
        "rna": "expression",
        "rnaseq": "expression",
        "rna_seq": "expression",
        "tpm": "expression",
    }
    return aliases.get(assay, assay)


def load_gencode_tss(path: Path) -> pd.DataFrame:
    rows = []
    with open_text(path) as handle:
        for line in handle:
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
            strand = parts[6]
            start_1 = int(parts[3])
            end_1 = int(parts[4])
            tss0 = start_1 - 1 if strand == "+" else end_1 - 1
            rows.append(
                {
                    "chr": chrom,
                    "tss0": tss0,
                    "strand": strand,
                    "gene_id": attrs.get("gene_id", "").split(".")[0],
                    "gene_name": attrs.get("gene_name", ""),
                    "source": "GENCODE",
                    "source_priority": 2,
                }
            )
    return pd.DataFrame(rows)


def load_fantom_tss_bed(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 6:
        raise ValueError("FANTOM TSS BED must have at least 6 columns: chr start end name score strand")
    rows = []
    for _, row in df.iterrows():
        chrom, start, end, name, _score, strand = row.iloc[:6]
        tss0 = int(start) if strand == "+" else int(end) - 1
        rows.append(
            {
                "chr": str(chrom) if str(chrom).startswith("chr") else f"chr{chrom}",
                "tss0": tss0,
                "strand": str(strand),
                "gene_id": str(name).split("|")[0],
                "gene_name": str(name).split("|")[-1],
                "source": "FANTOM5",
                "source_priority": 1,
            }
        )
    return pd.DataFrame(rows)


def choose_tss_sources(frames: list[pd.DataFrame]) -> pd.DataFrame:
    if not frames:
        raise ValueError("No TSS sources supplied.")
    tss = pd.concat(frames, ignore_index=True)
    tss = tss[tss["strand"].isin(["+", "-"])].copy()
    tss["gene_key"] = tss["gene_id"].fillna("").astype(str)
    tss.loc[tss["gene_key"] == "", "gene_key"] = tss["gene_name"].fillna("").astype(str)
    tss = tss.sort_values(["source_priority", "source", "chr", "tss0"])
    # Prefer FANTOM for duplicate TSS positions and for duplicate gene keys when present.
    tss = tss.drop_duplicates(["chr", "tss0", "strand"], keep="first")
    fantom_gene_keys = set(tss[(tss["source"] == "FANTOM5") & (tss["gene_key"] != "")]["gene_key"])
    keep = (tss["source"] == "FANTOM5") | (tss["gene_key"] == "") | (~tss["gene_key"].isin(fantom_gene_keys))
    return tss[keep].drop(columns=["source_priority"]).reset_index(drop=True)


def extract_sequence(fasta, chrom: str, start: int, end: int, strand: str) -> str | None:
    try:
        seq = str(fasta[chrom][start:end]).upper()
    except Exception:
        return None
    if strand == "-":
        seq = revcomp(seq)
    return seq if validate_seq(seq, end - start) else None


def build_promoter_windows(
    fasta_path: Path,
    gencode_gtf: Path | None,
    fantom_tss_bed: Path | None,
    upstream: int,
    downstream: int,
    max_promoters: int | None,
) -> pd.DataFrame:
    try:
        from pyfaidx import Fasta
    except ImportError as exc:
        raise SystemExit("Install pyfaidx first: python -m pip install pyfaidx") from exc

    frames = []
    if fantom_tss_bed:
        frames.append(load_fantom_tss_bed(fantom_tss_bed))
    if gencode_gtf:
        frames.append(load_gencode_tss(gencode_gtf))
    tss = choose_tss_sources(frames)
    if max_promoters:
        tss = tss.head(max_promoters)
    fasta = Fasta(str(fasta_path))
    rows = []
    for _, row in tss.iterrows():
        tss0 = int(row["tss0"])
        if row["strand"] == "+":
            start = max(0, tss0 - upstream)
            end = tss0 + downstream
        else:
            start = max(0, tss0 - downstream)
            end = tss0 + upstream
        seq = extract_sequence(fasta, row["chr"], start, end, row["strand"])
        if not seq:
            continue
        gene_key = str(row.get("gene_id") or row.get("gene_name") or f"{row['chr']}:{tss0}:{row['strand']}")
        rows.append(
            {
                "promoter_id": stable_promoter_id(row["chr"], start, end, row["strand"], gene_key),
                "chr": row["chr"],
                "start": start,
                "end": end,
                "strand": row["strand"],
                "tss0": tss0,
                "gene_id": row.get("gene_id", ""),
                "gene_name": row.get("gene_name", ""),
                "source": row.get("source", ""),
                "sequence": seq,
            }
        )
    return pd.DataFrame(rows).drop_duplicates("promoter_id")


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def validate_activity_rows(df: pd.DataFrame) -> pd.DataFrame:
    missing = REQUIRED_ACTIVITY_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"activity table missing required columns: {sorted(missing)}")
    out = df.copy()
    out["cell_type"] = out["cell_type"].map(harmonize_cell)
    out["assay"] = out["assay"].map(harmonize_assay)
    bad_cells = sorted(set(out["cell_type"]) - VALID_CELL_TYPES)
    bad_assays = sorted(set(out["assay"]) - VALID_ASSAYS)
    if bad_cells:
        raise ValueError(f"invalid cell_type values after harmonization: {bad_cells}")
    if bad_assays:
        raise ValueError(f"invalid assay values after harmonization: {bad_assays}")
    out["value"] = out["value"].map(normalize_activity_value)
    out = out.dropna(subset=["value"])
    return out


def bigwig_mean_signal(promoters: pd.DataFrame, manifest: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    try:
        import pyBigWig
    except ImportError as exc:
        raise SystemExit("Install pyBigWig first: python -m pip install pyBigWig") from exc

    man = read_tsv(manifest)
    required = {"cell_type", "assay", "path"}
    missing = required - set(man.columns)
    if missing:
        raise ValueError(f"signal manifest missing columns: {sorted(missing)}")
    rows = []
    source_rows = []
    for _, item in man.iterrows():
        path = Path(str(item["path"])).expanduser()
        if not path.exists():
            raise FileNotFoundError(path)
        cell = harmonize_cell(str(item["cell_type"]))
        assay = harmonize_assay(str(item["assay"]))
        if cell not in VALID_CELL_TYPES or assay not in VALID_ASSAYS:
            raise ValueError(f"invalid signal row cell/assay: {cell}/{assay}")
        bw = pyBigWig.open(str(path))
        try:
            for _, prom in promoters.iterrows():
                val = bw.stats(str(prom["chr"]), int(prom["start"]), int(prom["end"]), type="mean")[0]
                raw = 0.0 if val is None or math.isnan(float(val)) else float(val)
                rows.append(
                    {
                        "promoter_id": prom["promoter_id"],
                        "cell_type": cell,
                        "assay": assay,
                        "value": normalize_activity_value(raw),
                        "raw_value": raw,
                        "source": item.get("source", ""),
                        "accession": item.get("accession", ""),
                        "replicate": item.get("replicate", ""),
                    }
                )
        finally:
            bw.close()
        source_rows.append(
            {
                "cell_type": cell,
                "assay": assay,
                "source": item.get("source", ""),
                "accession": item.get("accession", ""),
                "replicate": item.get("replicate", ""),
                "path": str(path),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(source_rows)


def expression_activity(promoters: pd.DataFrame, expression_long: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = read_tsv(expression_long)
    required = {"cell_type", "value"}
    missing = required - set(expr.columns)
    if missing:
        raise ValueError(f"expression table missing columns: {sorted(missing)}")
    if "gene_id" not in expr.columns and "gene_name" not in expr.columns:
        raise ValueError("expression table must include gene_id or gene_name")
    expr = expr.copy()
    expr["cell_type"] = expr["cell_type"].map(harmonize_cell)
    expr["value"] = expr["value"].map(normalize_activity_value)
    expr = expr.dropna(subset=["value"])
    rows = []
    for key in ["gene_id", "gene_name"]:
        if key not in expr.columns or key not in promoters.columns:
            continue
        sub_expr = expr[expr[key].notna()].copy()
        sub_expr[key] = sub_expr[key].astype(str).str.split(".").str[0]
        prom = promoters[["promoter_id", key]].copy()
        prom[key] = prom[key].astype(str).str.split(".").str[0]
        merged = prom.merge(sub_expr, on=key, how="inner")
        for _, row in merged.iterrows():
            rows.append(
                {
                    "promoter_id": row["promoter_id"],
                    "cell_type": row["cell_type"],
                    "assay": "expression",
                    "value": row["value"],
                    "source": row.get("source", ""),
                    "accession": row.get("accession", ""),
                    "replicate": row.get("replicate", ""),
                }
            )
    source_cols = [c for c in ["cell_type", "source", "accession", "replicate"] if c in expr.columns]
    source_rows = expr[source_cols].drop_duplicates() if source_cols else pd.DataFrame()
    if not source_rows.empty:
        source_rows["assay"] = "expression"
        source_rows["path"] = str(expression_long)
    return pd.DataFrame(rows), source_rows


def tss_activity(promoters: pd.DataFrame, activity_path: Path, max_distance: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    tss = read_tsv(activity_path)
    if {"promoter_id", "cell_type", "value"} <= set(tss.columns):
        out = tss.copy()
        out["assay"] = "initiation"
        out = validate_activity_rows(out)
        return out, pd.DataFrame(
            [{"assay": "initiation", "source": "tss_activity_table", "path": str(activity_path)}]
        )
    required = {"chr", "tss0", "strand", "cell_type", "value"}
    missing = required - set(tss.columns)
    if missing:
        raise ValueError(f"TSS activity table missing columns: {sorted(missing)}")
    tss = tss.copy()
    tss["cell_type"] = tss["cell_type"].map(harmonize_cell)
    grouped: dict[tuple[str, str], pd.DataFrame] = {}
    positions: dict[tuple[str, str], list[int]] = {}
    for key, group in tss.groupby(["chr", "strand"]):
        group = group.sort_values("tss0").reset_index(drop=True)
        grouped[key] = group
        positions[key] = group["tss0"].astype(int).tolist()
    rows = []
    for _, prom in promoters.iterrows():
        key = (prom["chr"], prom["strand"])
        if key not in grouped:
            continue
        pos_list = positions[key]
        idx = bisect.bisect_left(pos_list, int(prom["tss0"]))
        candidate_idxs = [i for i in [idx - 1, idx, idx + 1] if 0 <= i < len(pos_list)]
        if not candidate_idxs:
            continue
        best = min(candidate_idxs, key=lambda i: abs(pos_list[i] - int(prom["tss0"])))
        if abs(pos_list[best] - int(prom["tss0"])) > max_distance:
            continue
        matches = grouped[key][grouped[key]["tss0"].astype(int) == pos_list[best]]
        for _, row in matches.iterrows():
            rows.append(
                {
                    "promoter_id": prom["promoter_id"],
                    "cell_type": row["cell_type"],
                    "assay": "initiation",
                    "value": normalize_activity_value(row["value"]),
                    "source": row.get("source", ""),
                    "accession": row.get("accession", ""),
                    "replicate": row.get("replicate", ""),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(
        [{"assay": "initiation", "source": "tss_activity_table", "path": str(activity_path)}]
    )


def median_combine_activity(frames: list[pd.DataFrame]) -> pd.DataFrame:
    frames = [f for f in frames if f is not None and not f.empty]
    if not frames:
        return pd.DataFrame(columns=sorted(REQUIRED_ACTIVITY_COLUMNS))
    df = pd.concat(frames, ignore_index=True, sort=False)
    df = validate_activity_rows(df)
    return (
        df.groupby(["promoter_id", "cell_type", "assay"], as_index=False)
        .agg(value=("value", "median"))
        .sort_values(["promoter_id", "cell_type", "assay"])
    )


def small_markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "No activity coverage."
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def write_qc_report(
    out: Path,
    promoters: pd.DataFrame,
    activity: pd.DataFrame,
    source_manifest: pd.DataFrame,
    strict: bool,
) -> dict[str, object]:
    coverage = (
        activity.groupby(["cell_type", "assay"])["promoter_id"].nunique().reset_index(name="promoters_with_signal")
        if not activity.empty
        else pd.DataFrame(columns=["cell_type", "assay", "promoters_with_signal"])
    )
    feature_grid = {(cell, assay) for cell in VALID_CELL_TYPES for assay in VALID_ASSAYS}
    observed = {(r["cell_type"], r["assay"]) for _, r in coverage.iterrows()}
    missing_grid = sorted([f"{cell}/{assay}" for cell, assay in feature_grid - observed])
    promoter_count = int(len(promoters))
    activity_count = int(len(activity))
    report = {
        "promoter_count": promoter_count,
        "activity_row_count": activity_count,
        "unique_activity_promoters": int(activity["promoter_id"].nunique()) if not activity.empty else 0,
        "cell_types": sorted(activity["cell_type"].unique().tolist()) if not activity.empty else [],
        "assays": sorted(activity["assay"].unique().tolist()) if not activity.empty else [],
        "missing_cell_assay_pairs": missing_grid,
        "source_count": int(len(source_manifest)),
        "strict": strict,
    }
    lines = [
        "# Stage 2 Real Data QC Report",
        "",
        f"- Promoters: {promoter_count:,}",
        f"- Activity rows: {activity_count:,}",
        f"- Promoters with any activity: {report['unique_activity_promoters']:,}",
        f"- Cell types: {', '.join(report['cell_types']) if report['cell_types'] else 'none'}",
        f"- Assays: {', '.join(report['assays']) if report['assays'] else 'none'}",
        f"- Source records: {len(source_manifest):,}",
        "",
        "## Coverage",
        "",
        small_markdown_table(coverage),
        "",
        "## Missing Cell/Assay Pairs",
        "",
        "\n".join(f"- {item}" for item in missing_grid) if missing_grid else "None.",
    ]
    out.write_text("\n".join(lines) + "\n")
    return report


def write_tsv(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def strict_coverage_check(activity: pd.DataFrame) -> None:
    required_cells = {"HSC", "HSPC", "T", "NK", "B", "MYELOID", "ERYTHROID", "MEGAKARYOCYTE"}
    observed_cells = set(activity["cell_type"]) if not activity.empty else set()
    missing_cells = sorted(required_cells - observed_cells)
    missing_assays = sorted(VALID_ASSAYS - set(activity["assay"])) if not activity.empty else sorted(VALID_ASSAYS)
    if missing_cells or missing_assays:
        raise SystemExit(f"Strict coverage failed. Missing cells={missing_cells}; missing assays={missing_assays}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=Path(__file__).resolve().parents[1] / "data" / "hsc_tnk_real")
    parser.add_argument("--reference-dir", default=Path(__file__).resolve().parents[1] / "downloads" / "references")
    parser.add_argument("--download-references", action="store_true")
    parser.add_argument("--hg38-fasta", default=None)
    parser.add_argument("--gencode-gtf", default=None)
    parser.add_argument("--fantom-tss-bed", default=None)
    parser.add_argument("--signal-manifest", default=None, help="TSV: cell_type, assay, path[,source,accession,replicate]")
    parser.add_argument("--expression-long", default=None, help="TSV: cell_type, gene_id/gene_name, value[,source,accession,replicate]")
    parser.add_argument("--tss-activity-long", default=None, help="TSV with promoter_id or chr/tss0/strand plus cell_type,value")
    parser.add_argument("--max-tss-distance", type=int, default=50)
    parser.add_argument("--upstream", type=int, default=150)
    parser.add_argument("--downstream", type=int, default=50)
    parser.add_argument("--max-promoters", type=int, default=None)
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--promoters-only", action="store_true")
    parser.add_argument("--skip-checksums", action="store_true")
    args = parser.parse_args()

    out_dir = ensure_dir(args.output_dir)
    reference_dir = ensure_dir(args.reference_dir)
    reports_dir = ensure_dir(Path(__file__).resolve().parents[1] / "reports" / "publishable")

    fasta_path = Path(args.hg38_fasta).expanduser() if args.hg38_fasta else None
    gtf_path = Path(args.gencode_gtf).expanduser() if args.gencode_gtf else None
    if args.download_references:
        hg38_gz = download(HG38_FASTA_URL, reference_dir / Path(urlparse(HG38_FASTA_URL).path).name)
        fasta_path = maybe_gunzip(hg38_gz)
        gtf_path = download(GENCODE_V38_URL, reference_dir / Path(urlparse(GENCODE_V38_URL).path).name)
    if not fasta_path or not fasta_path.exists():
        raise SystemExit("Provide --hg38-fasta or use --download-references.")
    if not gtf_path and not args.fantom_tss_bed:
        raise SystemExit("Provide --gencode-gtf, --fantom-tss-bed, or use --download-references.")

    fantom_path = Path(args.fantom_tss_bed).expanduser() if args.fantom_tss_bed else None
    promoters = build_promoter_windows(
        fasta_path=fasta_path,
        gencode_gtf=gtf_path,
        fantom_tss_bed=fantom_path,
        upstream=args.upstream,
        downstream=args.downstream,
        max_promoters=args.max_promoters,
    )
    write_tsv(out_dir / "promoter_windows.tsv", promoters)

    activity_frames: list[pd.DataFrame] = []
    harmonization_frames: list[pd.DataFrame] = []
    source_rows = [
        source_row(fasta_path, "reference_fasta", "UCSC hg38", skip_checksum=args.skip_checksums),
    ]
    if gtf_path:
        source_rows.append(source_row(gtf_path, "gene_annotation", "GENCODE v38", skip_checksum=args.skip_checksums))
    if fantom_path:
        source_rows.append(source_row(fantom_path, "tss_annotation", "FANTOM5", skip_checksum=args.skip_checksums))

    if args.signal_manifest:
        signal_df, signal_sources = bigwig_mean_signal(promoters, Path(args.signal_manifest).expanduser())
        activity_frames.append(signal_df)
        harmonization_frames.append(signal_sources)
        source_rows.append(source_row(Path(args.signal_manifest), "signal_manifest", skip_checksum=args.skip_checksums))
    if args.expression_long:
        expr_df, expr_sources = expression_activity(promoters, Path(args.expression_long).expanduser())
        activity_frames.append(expr_df)
        harmonization_frames.append(expr_sources)
        source_rows.append(source_row(Path(args.expression_long), "expression_long", skip_checksum=args.skip_checksums))
    if args.tss_activity_long:
        tss_df, tss_sources = tss_activity(promoters, Path(args.tss_activity_long).expanduser(), args.max_tss_distance)
        activity_frames.append(tss_df)
        harmonization_frames.append(tss_sources)
        source_rows.append(source_row(Path(args.tss_activity_long), "tss_activity_long", skip_checksum=args.skip_checksums))

    activity = median_combine_activity(activity_frames)
    if activity.empty and not args.promoters_only:
        raise SystemExit(
            "No activity rows were built. Provide --signal-manifest, --expression-long, "
            "--tss-activity-long, or pass --promoters-only."
        )
    if args.strict:
        strict_coverage_check(activity)
    write_tsv(out_dir / "activity_long.tsv", activity)

    source_manifest = pd.DataFrame(source_rows)
    write_tsv(out_dir / "source_manifest.lock.tsv", source_manifest)
    harmonization = pd.concat([f for f in harmonization_frames if f is not None and not f.empty], ignore_index=True, sort=False) if harmonization_frames else pd.DataFrame()
    write_tsv(out_dir / "sample_harmonization.tsv", harmonization)
    qc = write_qc_report(reports_dir / "data_qc_report.md", promoters, activity, source_manifest, args.strict)
    write_json(reports_dir / "data_qc_report.json", qc)

    print(f"Wrote real promoter windows to {out_dir / 'promoter_windows.tsv'}")
    print(f"Wrote real activity table to {out_dir / 'activity_long.tsv'}")
    print(f"Wrote QC report to {reports_dir / 'data_qc_report.md'}")


if __name__ == "__main__":
    main()
