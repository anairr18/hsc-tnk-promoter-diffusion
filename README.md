# HSC/T-NK Promoter Diffusion

Computational workflow for a two-stage promoter-generation project:

1. Finish a K562/HepG2/GM12878 DNA-Diffusion proof of concept.
2. Build the HSC-to-T/NK CAR promoter-design backbone using curated public hematopoietic inputs.

Wet-lab validation is not included. MPRA library design is included as a computational endpoint.

## Fastest Colab Path

First validate the software on a T4:

```python
import os, subprocess, tempfile
from pathlib import Path

BASE = Path(tempfile.gettempdir()) / "hsc_tnk_t4_smoke"
REPO = BASE / "hsc-tnk-promoter-diffusion"
BASE.mkdir(parents=True, exist_ok=True)
os.chdir(BASE)

if REPO.exists():
    subprocess.run(["git", "-C", str(REPO), "pull"], check=True)
else:
    subprocess.run(["git", "clone", "https://github.com/anairr18/hsc-tnk-promoter-diffusion.git", str(REPO)], check=True)

os.chdir(REPO)

env = os.environ.copy()
env["DNA_DIFFUSION_REPO"] = str(BASE / "DNA-Diffusion")
env["PATCH_PRETRAINED"] = "1"
env["SMOKE_EPOCHS"] = "3"
env["SMOKE_MIN_EPOCHS"] = "1"

subprocess.run(["bash", "run_t4_smoke.sh"], env=env, check=True)
print(Path("reports/t4_smoke_report.json").read_text())
```

Then run the full computational workflow on A100:

```python
import os, subprocess, tempfile
from pathlib import Path

BASE = Path(tempfile.gettempdir()) / "hsc_tnk_a100_full"
REPO = BASE / "hsc-tnk-promoter-diffusion"
BASE.mkdir(parents=True, exist_ok=True)
os.chdir(BASE)

if REPO.exists():
    subprocess.run(["git", "-C", str(REPO), "pull"], check=True)
else:
    subprocess.run(["git", "clone", "https://github.com/anairr18/hsc-tnk-promoter-diffusion.git", str(REPO)], check=True)

os.chdir(REPO)

env = os.environ.copy()
env["PATCH_PRETRAINED"] = "1"
env["DNA_DIFFUSION_REPO"] = str(BASE / "DNA-Diffusion")
env["DNA_DIFFUSION_CACHE"] = str(BASE / "dna_diffusion_cache")
env["DNA_DIFFUSION_EPOCHS"] = "5000"
env["DNA_DIFFUSION_MIN_EPOCHS"] = "2000"
env["STAGE2_MODE"] = "real"
env["BUILD_REAL_STAGE2_INPUTS"] = "1"
env["BUILD_ENCODE_STAGE2_DATA"] = "1"
env["BUILD_HCA_EXPRESSION"] = "1"
env["BUILD_HPA_EXPRESSION"] = "1"
env["DOWNLOAD_ENCODE_FILES"] = "1"
env["DOWNLOAD_REFERENCES"] = "1"
env["STAGE2_EPOCHS"] = "500"
env["STAGE2_MIN_EPOCHS"] = "100"
env["STAGE2_SAMPLES"] = "1000"

subprocess.run(["bash", "run_all_computational.sh"], env=env, check=True)
print(Path("reports/all_computational_summary.md").read_text())
```

For a cheap software-only check, set `STAGE2_MODE=demo`. Demo mode validates the plumbing but is not biological evidence.

## Curated Public Inputs

Run:

```bash
python scripts/curate_stage2_inputs_from_web.py --live-encode --limit-per-encode-query 25
```

Outputs:

- `reference_data/hematopoietic_manifest/curated_stage2_inputs.tsv`
- `reference_data/hematopoietic_manifest/curated_cell_state_source_map.tsv`
- `reference_data/hematopoietic_manifest/encode_live_experiment_candidates.tsv`
- `reference_data/hematopoietic_manifest/real_activity_long_schema.tsv`

The curated inputs cover FANTOM5 CAGE/TSS, ENCODE RAMPAGE/CAGE/RNA/accessibility, DICE immune TPM, BLUEPRINT hematopoietic epigenomes, Roadmap Epigenomics, HCA CD34 marrow single-cell expression, ENCODE SCREEN cCREs, JASPAR, HOCOMOCO, GENCODE v38, and hg38.

## Publishable Stage 2 Input Contract

For publishable HSC/T-NK candidate claims, build or provide real public-data inputs. The one-command public-data path is:

```bash
STAGE2_MODE=real \
BUILD_REAL_STAGE2_INPUTS=1 \
BUILD_ENCODE_STAGE2_DATA=1 \
DOWNLOAD_ENCODE_FILES=1 \
DOWNLOAD_REFERENCES=1 \
BUILD_HCA_EXPRESSION=1 \
BUILD_HPA_EXPRESSION=1 \
bash run_all_computational.sh
```

This automatically selects and downloads ENCODE GRCh38 bigWigs for accessibility, initiation, and promoter-proximal RNA signal, downloads a public marrow HCA/CELLxGENE H5AD for single-cell pseudobulk expression, and falls back to Human Protein Atlas single-cell type expression if the HCA asset is unavailable.

Optional extra public inputs can still be supplied:

```bash
FANTOM_TSS_BED=/path/to/fantom_hg38_tss.bed
SIGNAL_MANIFEST=/path/to/additional_signal_manifest.tsv
EXPRESSION_LONG_TSV=/path/to/additional_expression_long.tsv
TSS_ACTIVITY_LONG_TSV=/path/to/additional_tss_activity_long.tsv
```

If you already built the tables elsewhere, provide them directly:

```bash
STAGE2_MODE=real \
PROMOTER_WINDOWS=/path/to/promoter_windows.tsv \
ACTIVITY_LONG_TSV=/path/to/activity_long.tsv \
bash run_all_computational.sh
```

`PROMOTER_WINDOWS` columns:

```text
promoter_id chr start end strand tss0 gene_id gene_name source sequence
```

`ACTIVITY_LONG_TSV` columns:

```text
promoter_id cell_type assay value
```

Use harmonized `cell_type` values: `HSC`, `HSPC`, `T`, `NK`, `B`, `MYELOID`, `ERYTHROID`, `MEGAKARYOCYTE`.
Use assays: `accessibility`, `initiation`, `expression`.

`SIGNAL_MANIFEST` columns:

```text
cell_type assay path source accession replicate
```

Each `path` should point to an ATAC/DNase/CAGE/RAMPAGE bigWig. `assay` can be `accessibility` or `initiation`.

`EXPRESSION_LONG_TSV` columns:

```text
cell_type gene_id gene_name value source accession replicate
```

Use gene-level TPM-like values from DICE/HCA/ENCODE/BLUEPRINT; the builder applies `log1p`.

`TSS_ACTIVITY_LONG_TSV` columns:

```text
chr tss0 strand cell_type value source accession replicate
```

or:

```text
promoter_id cell_type value source accession replicate
```

Use FANTOM5 CAGE or ENCODE CAGE/RAMPAGE TSS activity. The builder applies `log1p` and joins nearest same-strand TSS within 50bp.

To preview ENCODE coverage without downloading large files:

```bash
python scripts/download_encode_stage2_public_data.py \
  --refresh-query \
  --output-dir downloads/encode_stage2_preview
```

To convert public HCA/BLUEPRINT AnnData into expression input:

```bash
python scripts/expression_from_anndata.py \
  --h5ad /path/to/public_bone_marrow_or_hematopoiesis.h5ad \
  --cell-type-column cell_type \
  --output data/hsc_tnk_real/hca_expression_long.tsv
```

To only inspect the automatically selected HCA/CELLxGENE asset without downloading the large H5AD:

```bash
python scripts/download_hca_stage2_expression.py --metadata-only
```

After a real run:

```bash
python scripts/validate_publishable_package.py --require-stage1
```

## Outputs

- Stage 1 report: `reports/cellline_poc_report.md`
- Stage 2 selected promoters: `outputs/hsc_tnk/tnk_specific_training_set.tsv`
- Stage 2 generated candidates: `outputs/hsc_tnk/generated/`
- Ranked candidates: `outputs/hsc_tnk/ranked_candidates.tsv`
- MPRA design table: `outputs/hsc_tnk/mpra_tnk_promoters.library.tsv`

See `IMPLEMENTATION_STATUS.md` and `JOURNAL_GRADE_BACKBONE.md` for the project-level notes.
