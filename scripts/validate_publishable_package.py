#!/usr/bin/env python3
"""Validate that outputs are suitable for computational manuscript claims."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from project_utils import validate_seq, write_json


REQUIRED_CELLS = {"HSC", "HSPC", "T", "NK", "B", "MYELOID", "ERYTHROID", "MEGAKARYOCYTE"}
REQUIRED_ASSAYS = {"accessibility", "initiation", "expression"}


def check_file(path: Path, failures: list[str], label: str) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        failures.append(f"Missing or empty {label}: {path}")
        return False
    return True


def no_demo_rows(path: Path, failures: list[str]) -> None:
    text = path.read_text(errors="replace").lower()
    if "demo" in text:
        failures.append(f"Demo marker found in publishable input/output: {path}")


def validate_activity(path: Path, failures: list[str], warnings: list[str]) -> dict[str, object]:
    df = pd.read_csv(path, sep="\t")
    required = {"promoter_id", "cell_type", "assay", "value"}
    missing = required - set(df.columns)
    if missing:
        failures.append(f"activity_long missing columns: {sorted(missing)}")
        return {}
    cells = set(df["cell_type"].astype(str))
    assays = set(df["assay"].astype(str))
    missing_cells = sorted(REQUIRED_CELLS - cells)
    missing_assays = sorted(REQUIRED_ASSAYS - assays)
    if missing_cells:
        failures.append(f"activity_long missing required cell states: {missing_cells}")
    if missing_assays:
        failures.append(f"activity_long missing required assays: {missing_assays}")
    if df["value"].isna().any():
        failures.append("activity_long contains NaN values")
    coverage = df.groupby(["cell_type", "assay"])["promoter_id"].nunique().reset_index(name="n_promoters")
    sparse = coverage[coverage["n_promoters"] < 100]
    if len(sparse):
        warnings.append("Some cell/assay pairs have fewer than 100 promoters with signal.")
    return {
        "activity_rows": int(len(df)),
        "activity_promoters": int(df["promoter_id"].nunique()),
        "cells": sorted(cells),
        "assays": sorted(assays),
    }


def validate_promoters(path: Path, failures: list[str]) -> dict[str, object]:
    df = pd.read_csv(path, sep="\t")
    required = {"promoter_id", "chr", "start", "end", "strand", "gene_id", "gene_name", "sequence"}
    missing = required - set(df.columns)
    if missing:
        failures.append(f"promoter_windows missing columns: {sorted(missing)}")
        return {}
    bad = int((~df["sequence"].map(validate_seq)).sum())
    if bad:
        failures.append(f"promoter_windows contains {bad} invalid non-200bp A/C/G/T sequences")
    if df["promoter_id"].duplicated().any():
        failures.append("promoter_windows contains duplicated promoter_id values")
    return {
        "promoters": int(len(df)),
        "valid_sequences": int(len(df) - bad),
        "sources": sorted(df.get("source", pd.Series(dtype=str)).dropna().astype(str).unique().tolist()),
    }


def validate_candidates(path: Path, failures: list[str], min_candidates: int) -> dict[str, object]:
    df = pd.read_csv(path, sep="\t")
    if "sequence" not in df.columns:
        failures.append("ranked_candidates lacks sequence column")
        return {}
    valid = df["sequence"].map(validate_seq)
    if int(valid.sum()) < min_candidates:
        failures.append(f"Only {int(valid.sum())} valid candidates; required >= {min_candidates}")
    if df["sequence"].duplicated().any():
        failures.append("ranked_candidates contains duplicated sequences")
    return {"ranked_candidates": int(len(df)), "valid_ranked_candidates": int(valid.sum())}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=Path(__file__).resolve().parents[1])
    parser.add_argument("--promoter-windows", default=None)
    parser.add_argument("--activity-long", default=None)
    parser.add_argument("--ranked-candidates", default=None)
    parser.add_argument("--mpra-library", default=None)
    parser.add_argument("--min-candidates", type=int, default=500)
    parser.add_argument("--require-stage1", action="store_true")
    parser.add_argument("--output-json", default=None)
    args = parser.parse_args()

    root = Path(args.project_root)
    failures: list[str] = []
    warnings: list[str] = []
    summary: dict[str, object] = {}

    promoter_path = Path(args.promoter_windows) if args.promoter_windows else root / "data" / "hsc_tnk_real" / "promoter_windows.tsv"
    activity_path = Path(args.activity_long) if args.activity_long else root / "data" / "hsc_tnk_real" / "activity_long.tsv"
    qc_report = root / "reports" / "publishable" / "data_qc_report.md"
    ranked_path = Path(args.ranked_candidates) if args.ranked_candidates else root / "outputs" / "hsc_tnk" / "ranked_candidates.tsv"
    library_path = Path(args.mpra_library) if args.mpra_library else root / "outputs" / "hsc_tnk" / "mpra_tnk_promoters.library.tsv"

    if check_file(promoter_path, failures, "real promoter windows"):
        no_demo_rows(promoter_path, failures)
        summary.update(validate_promoters(promoter_path, failures))
    if check_file(activity_path, failures, "real activity table"):
        no_demo_rows(activity_path, failures)
        summary.update(validate_activity(activity_path, failures, warnings))
    check_file(qc_report, failures, "data QC report")
    if check_file(ranked_path, failures, "ranked candidates"):
        no_demo_rows(ranked_path, failures)
        summary.update(validate_candidates(ranked_path, failures, args.min_candidates))
    check_file(library_path, failures, "MPRA library")

    if args.require_stage1:
        check_file(root / "reports" / "cellline_poc_report.md", failures, "Stage 1 POC report")

    result = {
        "pass": not failures,
        "failures": failures,
        "warnings": warnings,
        "summary": summary,
    }
    out_json = Path(args.output_json) if args.output_json else root / "reports" / "publishable" / "package_validation.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    write_json(out_json, result)
    md = out_json.with_suffix(".md")
    lines = [
        "# Publishable Package Validation",
        "",
        f"- Pass: {result['pass']}",
        "",
        "## Failures",
        "",
        "\n".join(f"- {item}" for item in failures) if failures else "None.",
        "",
        "## Warnings",
        "",
        "\n".join(f"- {item}" for item in warnings) if warnings else "None.",
    ]
    md.write_text("\n".join(lines) + "\n")
    print(json.dumps(result, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
