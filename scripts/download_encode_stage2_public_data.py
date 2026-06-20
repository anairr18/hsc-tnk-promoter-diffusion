#!/usr/bin/env python3
"""Download ENCODE public files for the real HSC/T-NK Stage 2 matrix.

Outputs a `signal_manifest.tsv` compatible with `build_real_stage2_inputs.py`.
The manifest can include accessibility, initiation, and expression bigWigs.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from urllib.parse import urlencode

import pandas as pd
import requests

from project_utils import ensure_dir, write_json


ENCODE_BASE = "https://www.encodeproject.org"

DEFAULT_QUERIES = [
    ("HSC", "CD34", "DNase-seq", ["cd34", "hematopoietic", "progenitor"]),
    ("HSC", "CD34", "CAGE", ["cd34", "hematopoietic", "progenitor"]),
    ("HSC", "CD34", "total RNA-seq", ["cd34", "hematopoietic", "progenitor"]),
    ("HSC", "CD34", "polyA plus RNA-seq", ["cd34", "hematopoietic", "progenitor"]),
    ("T", "T cell", "ATAC-seq", ["t cell", "t-helper", "cd4", "cd8", "thymus"]),
    ("T", "T cell", "DNase-seq", ["t cell", "t-helper", "cd4", "cd8", "thymus"]),
    ("T", "T cell", "total RNA-seq", ["t cell", "t-helper", "cd4", "cd8", "thymus"]),
    ("T", "T cell", "polyA plus RNA-seq", ["t cell", "t-helper", "cd4", "cd8", "thymus"]),
    ("NK", "natural killer", "ATAC-seq", ["natural killer", "cd56"]),
    ("NK", "natural killer", "DNase-seq", ["natural killer", "cd56"]),
    ("NK", "natural killer", "total RNA-seq", ["natural killer", "cd56"]),
    ("NK", "natural killer", "polyA plus RNA-seq", ["natural killer", "cd56"]),
    ("B", "B cell", "ATAC-seq", ["b cell", "b-cell", "lymphoblastoid"]),
    ("B", "B cell", "DNase-seq", ["b cell", "b-cell", "lymphoblastoid"]),
    ("B", "B cell", "total RNA-seq", ["b cell", "b-cell", "lymphoblastoid"]),
    ("B", "B cell", "polyA plus RNA-seq", ["b cell", "b-cell", "lymphoblastoid"]),
    ("MYELOID", "monocyte", "ATAC-seq", ["monocyte", "myeloid", "granulocyte"]),
    ("MYELOID", "monocyte", "DNase-seq", ["monocyte", "myeloid", "granulocyte"]),
    ("MYELOID", "monocyte", "total RNA-seq", ["monocyte", "myeloid", "granulocyte"]),
    ("MYELOID", "monocyte", "polyA plus RNA-seq", ["monocyte", "myeloid", "granulocyte"]),
    ("ERYTHROID", "erythroblast", "DNase-seq", ["erythroblast"]),
    ("ERYTHROID", "K562", "CAGE", ["k562"]),
    ("ERYTHROID", "K562", "total RNA-seq", ["k562"]),
    ("ERYTHROID", "K562", "polyA plus RNA-seq", ["k562"]),
]

ASSAY_TO_FEATURE = {
    "ATAC-seq": "accessibility",
    "DNase-seq": "accessibility",
    "CAGE": "initiation",
    "RAMPAGE": "initiation",
    "total RNA-seq": "expression",
    "polyA plus RNA-seq": "expression",
}

OUTPUT_PREFERENCE = {
    "accessibility": [
        "read-depth normalized signal",
        "fold change over control",
        "signal p-value",
    ],
    "initiation": [
        "plus strand signal of unique reads",
        "minus strand signal of unique reads",
        "plus strand signal of all reads",
        "minus strand signal of all reads",
    ],
    "expression": [
        "plus strand signal of unique reads",
        "minus strand signal of unique reads",
        "plus strand signal of all reads",
        "minus strand signal of all reads",
        "signal of unique reads",
        "signal of all reads",
    ],
}


def term_matches(term: str, text: str) -> bool:
    term = term.lower()
    text = text.lower()
    return bool(re.search(rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])", text))


def query_encode(cell_type: str, search_term: str, assay_title: str, include_terms: list[str], limit: int) -> list[dict]:
    params = {
        "type": "Experiment",
        "status": "released",
        "assay_title": assay_title,
        "searchTerm": search_term,
        "format": "json",
        "limit": str(limit),
    }
    url = f"{ENCODE_BASE}/search/?" + urlencode(params)
    data = requests.get(url, headers={"accept": "application/json"}, timeout=60).json()
    rows = []
    for experiment in data.get("@graph", []):
        biosample = experiment.get("biosample_ontology", {}).get("term_name", "")
        if include_terms and not any(term_matches(term, biosample) for term in include_terms):
            continue
        organisms = sorted(
            {
                rep.get("library", {}).get("biosample", {}).get("organism", {}).get("scientific_name", "")
                for rep in experiment.get("replicates", [])
            }
            - {""}
        )
        if organisms and "Homo sapiens" not in organisms:
            continue
        rows.append(
            {
                "cell_type": cell_type,
                "search_term": search_term,
                "assay_title": assay_title,
                "experiment_accession": experiment.get("accession", ""),
                "biosample": biosample,
                "organism": ",".join(organisms),
                "url": ENCODE_BASE + experiment.get("@id", ""),
            }
        )
    return rows


def fetch_experiment(accession: str) -> dict:
    response = requests.get(f"{ENCODE_BASE}/experiments/{accession}/?format=json", headers={"accept": "application/json"}, timeout=90)
    response.raise_for_status()
    return response.json()


def output_rank(feature: str, output_type: str, preferred_default: object) -> tuple[int, int]:
    prefs = OUTPUT_PREFERENCE[feature]
    try:
        output_idx = prefs.index(output_type)
    except ValueError:
        output_idx = 999
    preferred = 0 if preferred_default else 1
    return output_idx, preferred


def select_bigwigs(experiment: dict, feature: str, max_files: int) -> list[dict]:
    candidates = []
    for file_obj in experiment.get("files", []):
        if file_obj.get("status") != "released":
            continue
        if file_obj.get("file_format") != "bigWig":
            continue
        if file_obj.get("assembly") != "GRCh38":
            continue
        output_type = file_obj.get("output_type", "")
        if output_type not in OUTPUT_PREFERENCE[feature]:
            continue
        candidates.append(file_obj)
    candidates = sorted(
        candidates,
        key=lambda f: (
            output_rank(feature, f.get("output_type", ""), f.get("preferred_default")),
            int(f.get("file_size") or 0),
            f.get("accession", ""),
        ),
    )
    # Keep plus/minus pair for strand-specific assays when available; otherwise keep top files.
    if feature in {"initiation", "expression"}:
        plus = [f for f in candidates if "plus strand" in f.get("output_type", "")]
        minus = [f for f in candidates if "minus strand" in f.get("output_type", "")]
        selected = []
        if plus:
            selected.append(plus[0])
        if minus:
            selected.append(minus[0])
        if selected:
            return selected[:max_files]
    return candidates[:max_files]


def download_file(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        return
    with requests.get(url, stream=True, timeout=180) as response:
        response.raise_for_status()
        with path.open("wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                if chunk:
                    handle.write(chunk)


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=Path(__file__).resolve().parents[1] / "downloads" / "encode_stage2")
    parser.add_argument("--candidate-tsv", default=Path(__file__).resolve().parents[1] / "reference_data" / "hematopoietic_manifest" / "encode_live_experiment_candidates.tsv")
    parser.add_argument("--refresh-query", action="store_true", help="Query ENCODE live instead of using candidate TSV.")
    parser.add_argument("--download", action="store_true", help="Download selected bigWig files. Otherwise only writes URL manifests.")
    parser.add_argument("--limit-per-query", type=int, default=25)
    parser.add_argument("--max-experiments-per-cell-assay", type=int, default=2)
    parser.add_argument("--max-files-per-experiment", type=int, default=2)
    args = parser.parse_args()

    out = ensure_dir(args.output_dir)
    file_dir = ensure_dir(out / "files")
    meta_dir = ensure_dir(out / "metadata")

    if args.refresh_query:
        experiments = []
        for cell, term, assay, include in DEFAULT_QUERIES:
            experiments.extend(query_encode(cell, term, assay, include, args.limit_per_query))
    else:
        df = pd.read_csv(args.candidate_tsv, sep="\t")
        experiments = []
        for _, row in df.iterrows():
            assay = str(row.get("assay_title", ""))
            if assay not in ASSAY_TO_FEATURE:
                continue
            cell = str(row.get("harmonized_cell_state", row.get("cell_type", ""))).replace("HSC_HSPC", "HSC")
            if cell == "ERYTHROID" or cell == "MEGAKARYOCYTE":
                continue
            experiments.append(
                {
                    "cell_type": cell,
                    "search_term": row.get("search_term", ""),
                    "assay_title": assay,
                    "experiment_accession": row.get("accession", ""),
                    "biosample": row.get("biosample", ""),
                    "organism": row.get("organism", ""),
                    "url": row.get("url", ""),
                }
            )
        # Add broader lineage queries that the static curation may not include.
        for cell, term, assay, include in DEFAULT_QUERIES:
            if cell in {"ERYTHROID"}:
                experiments.extend(query_encode(cell, term, assay, include, args.limit_per_query))

    # De-duplicate experiments before file selection.
    seen = set()
    selected_experiments = []
    for row in experiments:
        acc = str(row.get("experiment_accession", ""))
        assay = str(row.get("assay_title", ""))
        if not acc or assay not in ASSAY_TO_FEATURE:
            continue
        feature = ASSAY_TO_FEATURE[assay]
        key = (row["cell_type"], feature)
        if (key, acc) in seen:
            continue
        seen.add((key, acc))
        selected_experiments.append(row)

    signal_rows = []
    selected_file_rows = []
    failed = []
    successful_experiments_by_key: dict[tuple[str, str], int] = {}
    for row in selected_experiments:
        acc = row["experiment_accession"]
        feature = ASSAY_TO_FEATURE[row["assay_title"]]
        key = (row["cell_type"], feature)
        if successful_experiments_by_key.get(key, 0) >= args.max_experiments_per_cell_assay:
            continue
        try:
            exp = fetch_experiment(acc)
        except Exception as exc:
            failed.append({**row, "error": str(exc)})
            continue
        (meta_dir / f"{acc}.json").write_text(json.dumps(exp, indent=2) + "\n")
        files = select_bigwigs(exp, feature, args.max_files_per_experiment)
        if not files:
            failed.append({**row, "error": "no usable GRCh38 bigWig selected"})
            continue
        successful_experiments_by_key[key] = successful_experiments_by_key.get(key, 0) + 1
        for file_obj in files:
            href = file_obj.get("href", "")
            url = ENCODE_BASE + href
            filename = Path(href).name
            local_path = file_dir / filename
            if args.download:
                download_file(url, local_path)
                manifest_path = str(local_path)
            else:
                manifest_path = url
            out_row = {
                "cell_type": row["cell_type"],
                "assay": feature,
                "path": manifest_path,
                "source": "ENCODE",
                "accession": acc,
                "replicate": file_obj.get("accession", ""),
                "file_accession": file_obj.get("accession", ""),
                "output_type": file_obj.get("output_type", ""),
                "assembly": file_obj.get("assembly", ""),
                "biosample": row.get("biosample", ""),
                "download_url": url,
                "bytes": file_obj.get("file_size", ""),
            }
            signal_rows.append(out_row)
            selected_file_rows.append(out_row)

    write_tsv(
        out / "signal_manifest.tsv",
        signal_rows,
        ["cell_type", "assay", "path", "source", "accession", "replicate", "file_accession", "output_type", "assembly", "biosample", "download_url", "bytes"],
    )
    write_tsv(
        out / "selected_encode_files.tsv",
        selected_file_rows,
        ["cell_type", "assay", "path", "source", "accession", "replicate", "file_accession", "output_type", "assembly", "biosample", "download_url", "bytes"],
    )
    write_tsv(
        out / "selected_encode_experiments.tsv",
        selected_experiments,
        ["cell_type", "search_term", "assay_title", "experiment_accession", "biosample", "organism", "url"],
    )
    write_tsv(out / "failed_encode_experiments.tsv", failed, ["cell_type", "search_term", "assay_title", "experiment_accession", "biosample", "organism", "url", "error"])
    coverage = pd.DataFrame(signal_rows).groupby(["cell_type", "assay"]).size().reset_index(name="files") if signal_rows else pd.DataFrame(columns=["cell_type", "assay", "files"])
    coverage.to_csv(out / "encode_signal_coverage.tsv", sep="\t", index=False)
    write_json(
        out / "encode_stage2_download_summary.json",
        {
            "downloaded_files": bool(args.download),
            "signal_manifest": str(out / "signal_manifest.tsv"),
            "selected_file_count": len(signal_rows),
            "selected_experiment_count": len(selected_experiments),
            "coverage": coverage.to_dict(orient="records"),
        },
    )
    print(f"Wrote ENCODE Stage 2 manifests to {out}")
    if not args.download:
        print("Use --download to fetch selected bigWigs; current manifest paths are URLs.")


if __name__ == "__main__":
    main()
