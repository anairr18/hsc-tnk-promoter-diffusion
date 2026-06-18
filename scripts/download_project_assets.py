#!/usr/bin/env python3
"""Download public assets and metadata needed by the project.

Large downloads are opt-in. By default this fetches lightweight metadata and
records exact URLs for large genome/omics assets so the run is reproducible.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import requests

from project_utils import ensure_dir


CELL_LINE_EXPERIMENTS = {
    "K562_DNase": "ENCSR000EOT",
    "K562_RNA_1": "ENCSR000AEQ",
    "K562_RNA_2": "ENCSR000COK",
    "HepG2_DNase": "ENCSR000ENQ",
    "HepG2_RNA_1": "ENCSR000EYR",
    "HepG2_RNA_2": "ENCSR931WGT",
    "GM12878_DNase": "ENCSR000EMT",
    "GM12878_RNA_1": "ENCSR000AEF",
    "GM12878_RNA_2": "ENCSR000AEG",
}

URLS = {
    "gencode_v38_gtf": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz",
    "hg38_fasta": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    "fantom5": "https://fantom.gsc.riken.jp/5/",
    "encode": "https://www.encodeproject.org/",
    "blueprint": "https://projects.ensembl.org/blueprint/",
    "roadmap": "https://egg2.wustl.edu/roadmap/web_portal/processed_data.html",
    "dice": "https://dice-database.org/",
    "hca": "https://data.humancellatlas.org/",
    "encode_screen": "https://screen.encodeproject.org/",
    "jaspar": "https://jaspar.elixir.no/",
}


def download(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with path.open("wb") as f:
            for chunk in r.iter_content(1024 * 1024):
                if chunk:
                    f.write(chunk)


def fetch_encode_metadata(out: Path) -> None:
    meta_dir = ensure_dir(out / "encode_metadata")
    for label, accession in CELL_LINE_EXPERIMENTS.items():
        url = f"https://www.encodeproject.org/experiments/{accession}/?format=json"
        r = requests.get(url, headers={"accept": "application/json"}, timeout=60)
        r.raise_for_status()
        (meta_dir / f"{label}_{accession}.json").write_text(json.dumps(r.json(), indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default=Path(__file__).resolve().parents[1] / "downloads")
    parser.add_argument("--download-gencode", action="store_true")
    parser.add_argument("--download-hg38", action="store_true", help="Large: about 1 GB compressed.")
    parser.add_argument("--skip-encode-metadata", action="store_true")
    args = parser.parse_args()

    out = ensure_dir(args.output_dir)
    (out / "asset_urls.json").write_text(json.dumps(URLS, indent=2) + "\n")
    if not args.skip_encode_metadata:
        fetch_encode_metadata(out)
    if args.download_gencode:
        download(URLS["gencode_v38_gtf"], out / "gencode.v38.annotation.gtf.gz")
    if args.download_hg38:
        download(URLS["hg38_fasta"], out / "hg38.fa.gz")
    print(f"Wrote asset manifest and requested downloads to {out}")


if __name__ == "__main__":
    main()
