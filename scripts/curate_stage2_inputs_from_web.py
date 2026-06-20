#!/usr/bin/env python3
"""Curate public inputs for the HSC/T-NK promoter-design stage.

The static manifest is based on public source pages for FANTOM5, ENCODE,
BLUEPRINT, DICE, HCA, Roadmap, SCREEN, JASPAR, and HOCOMOCO. The optional
ENCODE query adds live experiment accessions for the major blood lineages.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Any
from urllib.parse import urlencode

import requests

from project_utils import ensure_dir, write_json


CURATED_SOURCES: list[dict[str, str]] = [
    {
        "priority": "1",
        "source_id": "fantom5_grch38_cage_tss",
        "source": "FANTOM5 / UCSC",
        "url": "https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=fantom5",
        "modality": "CAGE TSS positions and sample TPM/counts",
        "assembly": "GRCh38/hg38",
        "access": "public",
        "cell_states": "human primary cells, immune cells, blood cells, cell lines",
        "stage2_use": "Preferred TSS source and transcription-initiation signal for promoter windows.",
        "download_action": "Download hg38 FANTOM5 DPI/TSS peaks and TSS activity TPM/count tracks.",
        "notes": "Use promoter/TSS activity as the initiation feature; map sample names into the harmonized lineage labels.",
    },
    {
        "priority": "2",
        "source_id": "encode_cage_rampage",
        "source": "ENCODE",
        "url": "https://www.encodeproject.org/data-standards/rampage/",
        "modality": "RAMPAGE/CAGE",
        "assembly": "mostly GRCh38/hg38",
        "access": "public",
        "cell_states": "available ENCODE biosamples, including selected hematopoietic samples",
        "stage2_use": "Supplement FANTOM5 with direct TSS/initiation evidence where matched cell states exist.",
        "download_action": "Query ENCODE for released RAMPAGE/CAGE bigWig or peak outputs.",
        "notes": "RAMPAGE is base-pair TSS evidence and quantifies TSS expression; use when cell-state matched.",
    },
    {
        "priority": "3",
        "source_id": "encode_accessibility_rnaseq",
        "source": "ENCODE",
        "url": "https://www.encodeproject.org/",
        "modality": "DNase-seq, ATAC-seq, total/polyA RNA-seq",
        "assembly": "mostly GRCh38/hg38",
        "access": "public",
        "cell_states": "CD34/HSPC, T, NK, B, monocyte/myeloid, erythroid, megakaryocyte when available",
        "stage2_use": "Accessibility and expression features for target/off-target promoter scoring.",
        "download_action": "Use live ENCODE query output from this script, then fetch preferred bigWig/TSV outputs.",
        "notes": "Prefer matched assay/cell-state data; otherwise use BLUEPRINT/Roadmap/DICE fallback.",
    },
    {
        "priority": "4",
        "source_id": "dice_immune_tpm",
        "source": "DICE",
        "url": "https://dice-database.org/downloads",
        "modality": "bulk RNA-seq TPM and immune eQTLs",
        "assembly": "GRCh37.p19 for blood dataset; gene-level expression",
        "access": "public",
        "cell_states": "naive B, monocytes, NK CD56dim CD16+, CD4/CD8/Treg/Th subsets",
        "stage2_use": "Mature immune expression reference for T/NK targets and B/myeloid off-targets.",
        "download_action": "Download all-cell-type mean TPM plus per-cell TPM tables.",
        "notes": "Gene-level expression should be joined to promoter genes; no coordinate liftover needed for TPM.",
    },
    {
        "priority": "5",
        "source_id": "blueprint_processed_hematopoietic",
        "source": "BLUEPRINT",
        "url": "https://projects.ensembl.org/blueprint/",
        "modality": "RNA-seq, DNase-seq, ChIP-seq/histone marks, methylation",
        "assembly": "source-dependent; processed files include reference metadata",
        "access": "mixed public/controlled",
        "cell_states": "hematopoietic lineages including monocyte examples and broader blood cell releases",
        "stage2_use": "Primary hematopoietic epigenome/expression backbone when accessible.",
        "download_action": "Use processed FTP files when public; record controlled-access datasets separately.",
        "notes": "High-value source for blood-cell epigenomes; not all samples are immediately downloadable.",
    },
    {
        "priority": "6",
        "source_id": "roadmap_epigenomics",
        "source": "NIH Roadmap Epigenomics",
        "url": "https://registry.opendata.aws/roadmapepigenomics/",
        "modality": "DNase, histone marks, methylation, mRNA expression",
        "assembly": "mostly hg19",
        "access": "public",
        "cell_states": "CD34 and immune/blood reference states where available",
        "stage2_use": "Fallback epigenome tracks and sanity-check promoter chromatin states.",
        "download_action": "Use AWS/Open Data directory listing; lift hg19 coordinates to hg38 if coordinate-level.",
        "notes": "Public and stable, but older assembly and older assays.",
    },
    {
        "priority": "7",
        "source_id": "hca_cd34_bone_marrow",
        "source": "Human Cell Atlas",
        "url": "https://explore.data.humancellatlas.org/projects/091cf39b-01bc-42e5-9437-f419a66c8a45",
        "modality": "single-cell RNA-seq",
        "assembly": "gene-level processed matrices",
        "access": "public",
        "cell_states": "CD34-positive, CD38-negative HSCs and bone marrow differentiation context",
        "stage2_use": "HSC/HSPC leakiness baseline and pseudobulk differentiation-state expression.",
        "download_action": "Download h5ad/loom processed outputs and build HSC/progenitor pseudobulk tables.",
        "notes": "Use as biological context and HSC expression baseline; not a direct promoter initiation assay.",
    },
    {
        "priority": "8",
        "source_id": "encode_screen_ccres",
        "source": "ENCODE SCREEN / UCSC",
        "url": "https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=encodeCcreCombined",
        "modality": "cCRE annotations and biosample-specific regulatory support",
        "assembly": "GRCh38/hg38",
        "access": "public",
        "cell_states": "ENCODE4 core biosamples and integrated cCRE registry",
        "stage2_use": "Promoter/enhancer annotation, active promoter support, and genomic QC.",
        "download_action": "Download integrated hg38 cCREs and biosample-specific cCRE support if needed.",
        "notes": "Use to tag promoters as PLS/pELS/dELS/CTCF-like and deprioritize ambiguous regions.",
    },
    {
        "priority": "9",
        "source_id": "jaspar_core_vertebrates",
        "source": "JASPAR",
        "url": "https://jaspar.elixir.no/downloads/",
        "modality": "TF motif PFMs/PWMs",
        "assembly": "not applicable",
        "access": "public",
        "cell_states": "not cell-state-specific",
        "stage2_use": "Motif enrichment, interpretability, target lineage motif checks, unwanted motif filters.",
        "download_action": "Download JASPAR CORE vertebrate non-redundant PFMs in MEME or JASPAR format.",
        "notes": "Use alongside HOCOMOCO; do not overfit ranking only to motif presence.",
    },
    {
        "priority": "10",
        "source_id": "hocomoco_v12_human",
        "source": "HOCOMOCO",
        "url": "https://hocomoco14.autosome.org/downloads_v12",
        "modality": "human/mouse TF motif PWMs/PFMs",
        "assembly": "not applicable",
        "access": "public",
        "cell_states": "not cell-state-specific",
        "stage2_use": "Secondary motif library for robustness and unwanted off-target motif filters.",
        "download_action": "Download H12CORE or H12INVIVO human PWM/PFM flat files.",
        "notes": "Use human subset for scanning 200bp promoter candidates.",
    },
    {
        "priority": "11",
        "source_id": "gencode_v38",
        "source": "GENCODE",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz",
        "modality": "gene annotation",
        "assembly": "GRCh38/hg38",
        "access": "public",
        "cell_states": "all genes",
        "stage2_use": "Fallback TSS and gene ID mapping when FANTOM5 TSS is absent.",
        "download_action": "Download gencode.v38.annotation.gtf.gz.",
        "notes": "Use protein-coding genes by default unless promoter atlas says otherwise.",
    },
    {
        "priority": "12",
        "source_id": "ucsc_hg38_fasta",
        "source": "UCSC",
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        "modality": "reference genome FASTA",
        "assembly": "GRCh38/hg38",
        "access": "public",
        "cell_states": "not applicable",
        "stage2_use": "Extract 200bp promoter-proximal windows.",
        "download_action": "Download and index hg38.fa.gz or use a mounted reference FASTA.",
        "notes": "Needed only when building promoter windows from coordinates.",
    },
]


CELL_STATE_SOURCE_MAP: list[dict[str, str]] = [
    {
        "harmonized_cell_state": "HSC",
        "role": "off-target leakiness baseline",
        "primary_sources": "HCA CD34 bone marrow; FANTOM5 CAGE if matched; BLUEPRINT/Roadmap CD34 tracks",
        "fallback_sources": "ENCODE CD34 DNase/RNA/CAGE candidates from live query",
        "required_assays": "expression, accessibility, initiation if available",
    },
    {
        "harmonized_cell_state": "HSPC",
        "role": "off-target leakiness baseline and progenitor specificity",
        "primary_sources": "HCA pseudobulk; BLUEPRINT/Roadmap progenitors; ENCODE CD34/progenitor assays",
        "fallback_sources": "GENCODE-linked expression from HCA only",
        "required_assays": "expression, accessibility",
    },
    {
        "harmonized_cell_state": "T",
        "role": "target activity",
        "primary_sources": "DICE T-cell TPM; FANTOM5 T-cell CAGE; ENCODE T-cell ATAC/DNase/RNA",
        "fallback_sources": "BLUEPRINT/Roadmap T-cell tracks",
        "required_assays": "expression, accessibility, initiation if available",
    },
    {
        "harmonized_cell_state": "NK",
        "role": "target activity",
        "primary_sources": "DICE NK TPM; ENCODE NK ATAC/DNase/RNA; FANTOM5 NK-like CAGE if available",
        "fallback_sources": "BLUEPRINT/Roadmap NK tracks",
        "required_assays": "expression, accessibility",
    },
    {
        "harmonized_cell_state": "B",
        "role": "off-target mature lymphoid",
        "primary_sources": "DICE B-cell TPM; ENCODE/Roadmap/BLUEPRINT B-cell tracks",
        "fallback_sources": "FANTOM5 B-cell CAGE",
        "required_assays": "expression, accessibility",
    },
    {
        "harmonized_cell_state": "MYELOID",
        "role": "off-target myeloid",
        "primary_sources": "DICE monocyte TPM; BLUEPRINT monocyte DNase/RNA; ENCODE monocyte assays",
        "fallback_sources": "Roadmap monocyte/macrophage tracks",
        "required_assays": "expression, accessibility",
    },
    {
        "harmonized_cell_state": "ERYTHROID",
        "role": "off-target erythroid lineage",
        "primary_sources": "ENCODE erythroid progenitor RNA/ATAC/DNase; BLUEPRINT/Roadmap erythroid tracks",
        "fallback_sources": "HCA marrow pseudobulk",
        "required_assays": "expression, accessibility",
    },
    {
        "harmonized_cell_state": "MEGAKARYOCYTE",
        "role": "off-target megakaryocyte lineage",
        "primary_sources": "ENCODE megakaryocyte RNA/ATAC; BLUEPRINT/Roadmap megakaryocyte tracks",
        "fallback_sources": "HCA marrow pseudobulk",
        "required_assays": "expression, accessibility",
    },
]


ENCODE_QUERIES: list[dict[str, Any]] = [
    {"harmonized_cell_state": "HSC_HSPC", "search_term": "CD34", "assay_title": "DNase-seq", "include_terms": ["cd34", "hematopoietic", "progenitor"]},
    {"harmonized_cell_state": "HSC_HSPC", "search_term": "CD34", "assay_title": "polyA plus RNA-seq", "include_terms": ["cd34", "hematopoietic", "progenitor"]},
    {"harmonized_cell_state": "HSC_HSPC", "search_term": "CD34", "assay_title": "total RNA-seq", "include_terms": ["cd34", "hematopoietic", "progenitor"]},
    {"harmonized_cell_state": "HSC_HSPC", "search_term": "CD34", "assay_title": "CAGE", "include_terms": ["cd34", "hematopoietic", "progenitor"]},
    {"harmonized_cell_state": "T", "search_term": "T cell", "assay_title": "ATAC-seq", "include_terms": ["t cell", "t-helper", "cd4", "cd8", "thymus"]},
    {"harmonized_cell_state": "T", "search_term": "T cell", "assay_title": "DNase-seq", "include_terms": ["t cell", "t-helper", "cd4", "cd8", "thymus"]},
    {"harmonized_cell_state": "T", "search_term": "T cell", "assay_title": "polyA plus RNA-seq", "include_terms": ["t cell", "t-helper", "cd4", "cd8", "thymus"]},
    {"harmonized_cell_state": "T", "search_term": "T cell", "assay_title": "total RNA-seq", "include_terms": ["t cell", "t-helper", "cd4", "cd8", "thymus"]},
    {"harmonized_cell_state": "NK", "search_term": "natural killer", "assay_title": "ATAC-seq", "include_terms": ["natural killer", "cd56"]},
    {"harmonized_cell_state": "NK", "search_term": "natural killer", "assay_title": "DNase-seq", "include_terms": ["natural killer", "cd56"]},
    {"harmonized_cell_state": "NK", "search_term": "natural killer", "assay_title": "polyA plus RNA-seq", "include_terms": ["natural killer", "cd56"]},
    {"harmonized_cell_state": "NK", "search_term": "natural killer", "assay_title": "total RNA-seq", "include_terms": ["natural killer", "cd56"]},
    {"harmonized_cell_state": "B", "search_term": "B cell", "assay_title": "ATAC-seq", "include_terms": ["b cell", "b-cell", "lymphoblastoid"]},
    {"harmonized_cell_state": "B", "search_term": "B cell", "assay_title": "DNase-seq", "include_terms": ["b cell", "b-cell", "lymphoblastoid"]},
    {"harmonized_cell_state": "B", "search_term": "B cell", "assay_title": "polyA plus RNA-seq", "include_terms": ["b cell", "b-cell", "lymphoblastoid"]},
    {"harmonized_cell_state": "MYELOID", "search_term": "monocyte", "assay_title": "ATAC-seq", "include_terms": ["monocyte", "myeloid", "granulocyte"]},
    {"harmonized_cell_state": "MYELOID", "search_term": "monocyte", "assay_title": "DNase-seq", "include_terms": ["monocyte", "myeloid", "granulocyte"]},
    {"harmonized_cell_state": "MYELOID", "search_term": "monocyte", "assay_title": "polyA plus RNA-seq", "include_terms": ["monocyte", "myeloid", "granulocyte"]},
    {"harmonized_cell_state": "ERYTHROID", "search_term": "erythroid", "assay_title": "ATAC-seq", "include_terms": ["erythroid", "megakaryocyte-erythroid"]},
    {"harmonized_cell_state": "ERYTHROID", "search_term": "erythroid", "assay_title": "total RNA-seq", "include_terms": ["erythroid", "megakaryocyte-erythroid"]},
    {"harmonized_cell_state": "MEGAKARYOCYTE", "search_term": "megakaryocyte", "assay_title": "ATAC-seq", "include_terms": ["megakaryocyte"]},
    {"harmonized_cell_state": "MEGAKARYOCYTE", "search_term": "megakaryocyte", "assay_title": "total RNA-seq", "include_terms": ["megakaryocyte"]},
]


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    if not rows:
        if fieldnames is None:
            fieldnames = []
        with path.open("w", newline="") as f:
            csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t").writeheader()
        return
    fieldnames = fieldnames or list(rows[0].keys())
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def term_matches_biosample(term: str, biosample: str) -> bool:
    term = term.lower()
    biosample = biosample.lower()
    if re.search(rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])", biosample):
        return True
    return False


def query_encode(limit_per_query: int = 10) -> list[dict[str, str]]:
    """Return lightweight released ENCODE experiment candidates."""
    headers = {"accept": "application/json"}
    rows: list[dict[str, str]] = []
    for query in ENCODE_QUERIES:
        params = {
            "type": "Experiment",
            "status": "released",
            "assay_title": query["assay_title"],
            "searchTerm": query["search_term"],
            "format": "json",
            "limit": str(limit_per_query),
        }
        url = "https://www.encodeproject.org/search/?" + urlencode(params)
        try:
            response = requests.get(url, headers=headers, timeout=45)
            data = response.json()
        except Exception as exc:  # pragma: no cover - network dependent
            public_query = {k: v for k, v in query.items() if k != "include_terms"}
            rows.append(
                {
                    **public_query,
                    "accession": "",
                    "biosample": "",
                    "status": "query_failed",
                    "total_for_query": "",
                    "url": url,
                    "notes": str(exc),
                }
            )
            continue
        for experiment in data.get("@graph", []):
            biosample = experiment.get("biosample_ontology", {}).get("term_name", "")
            include_terms = [str(term).lower() for term in query.get("include_terms", [])]
            if include_terms and not any(term_matches_biosample(term, biosample) for term in include_terms):
                continue
            organisms = sorted(
                {
                    rep.get("library", {})
                    .get("biosample", {})
                    .get("organism", {})
                    .get("scientific_name", "")
                    for rep in experiment.get("replicates", [])
                }
                - {""}
            )
            if organisms and "Homo sapiens" not in organisms:
                continue
            public_query = {k: v for k, v in query.items() if k != "include_terms"}
            rows.append(
                {
                    **public_query,
                    "accession": experiment.get("accession", ""),
                    "biosample": biosample,
                    "organism": ",".join(organisms) if organisms else "",
                    "status": experiment.get("status", ""),
                    "total_for_query": str(data.get("total", "")),
                    "url": "https://www.encodeproject.org" + experiment.get("@id", ""),
                    "notes": experiment.get("description", "") or "",
                }
            )
    deduped: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        key = (row.get("harmonized_cell_state", ""), row.get("assay_title", ""), row.get("accession", ""))
        deduped[key] = row
    return list(deduped.values())


def write_activity_schema(out: Path) -> None:
    rows = [
        {
            "promoter_id": "prom_00000001",
            "cell_type": "T",
            "assay": "accessibility",
            "value": "2.71",
            "recommended_source": "ENCODE/DICE/BLUEPRINT/Roadmap after promoter overlap or gene join",
            "notes": "Numeric, non-negative. Use log1p-normalized signal or TPM-like scale consistently.",
        },
        {
            "promoter_id": "prom_00000001",
            "cell_type": "NK",
            "assay": "initiation",
            "value": "4.13",
            "recommended_source": "FANTOM5 CAGE or ENCODE RAMPAGE/CAGE",
            "notes": "Promoter/TSS signal should be assigned to the 200bp window.",
        },
        {
            "promoter_id": "prom_00000001",
            "cell_type": "HSC",
            "assay": "expression",
            "value": "0.12",
            "recommended_source": "HCA/DICE/ENCODE/BLUEPRINT gene-level expression",
            "notes": "Join by gene_id/gene_name from promoter window table.",
        },
    ]
    write_tsv(out, rows)


def write_readme(out: Path, live_encode: bool) -> None:
    text = f"""# Curated Stage 2 Inputs

Generated by `scripts/curate_stage2_inputs_from_web.py`.

## Files

- `curated_stage2_inputs.tsv`: source-level manifest for real HSC/T-NK promoter design inputs.
- `curated_cell_state_source_map.tsv`: harmonized lineage labels and preferred sources.
- `real_activity_long_schema.tsv`: required schema for `scripts/merge_promoter_activity_matrix.py`.
- `encode_live_experiment_candidates.tsv`: optional live ENCODE experiment candidates. Present: `{str(live_encode)}`.

## Minimal Real Stage 2 Input Contract

To run `STAGE2_MODE=real`, provide:

- `PROMOTER_WINDOWS`: TSV with `promoter_id`, `chr`, `start`, `end`, `strand`, `gene_id`, `gene_name`, `sequence`.
- `ACTIVITY_LONG_TSV`: TSV with `promoter_id`, `cell_type`, `assay`, `value`.

Use harmonized `cell_type` values: `HSC`, `HSPC`, `T`, `NK`, `B`, `MYELOID`, `ERYTHROID`, `MEGAKARYOCYTE`.
Use assays: `accessibility`, `initiation`, `expression`.

The software can run in demo mode immediately. Real biological claims require replacing demo inputs with these curated public sources.
"""
    (out / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out-dir",
        default=Path(__file__).resolve().parents[1] / "reference_data" / "hematopoietic_manifest",
    )
    parser.add_argument("--live-encode", action="store_true", help="Query ENCODE for released experiment accessions.")
    parser.add_argument("--limit-per-encode-query", type=int, default=10)
    args = parser.parse_args()

    out = ensure_dir(args.out_dir)
    write_tsv(out / "curated_stage2_inputs.tsv", CURATED_SOURCES)
    write_json(out / "curated_stage2_inputs.json", CURATED_SOURCES)
    write_tsv(out / "curated_cell_state_source_map.tsv", CELL_STATE_SOURCE_MAP)
    write_json(out / "curated_cell_state_source_map.json", CELL_STATE_SOURCE_MAP)
    write_activity_schema(out / "real_activity_long_schema.tsv")

    if args.live_encode:
        rows = query_encode(args.limit_per_encode_query)
        write_tsv(
            out / "encode_live_experiment_candidates.tsv",
            rows,
            [
                "harmonized_cell_state",
                "search_term",
                "assay_title",
                "accession",
                "biosample",
                "organism",
                "status",
                "total_for_query",
                "url",
                "notes",
            ],
        )

    write_readme(out, args.live_encode)
    print(f"Wrote curated Stage 2 input manifests to {out}")


if __name__ == "__main__":
    main()
