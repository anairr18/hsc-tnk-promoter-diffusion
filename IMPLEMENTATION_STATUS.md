# Implementation Status

This export now has two runnable layers:

1. **Cell-line proof of concept** for K562, HepG2, and GM12878.
2. **HSC/T-NK CAR promoter pipeline scaffolding** for public-data training and BCH donor-cell MPRA validation.

## Stage 1: Cell-Line POC

Run dataset/training on a GPU machine:

```bash
python -m pip install -r requirements.txt
python pipeline/cell1_dataset_and_training.py
python pipeline/cell2_generate_and_backup.py
```

Then validate generated sequences against pretrained and endogenous controls:

```bash
python scripts/validate_cellline_poc.py --make-plots
```

For a stronger computational-journal POC, also create leakage-resistant split bundles and audit novelty:

```bash
python scripts/prepare_high_impact_poc.py \
  --input ~/dna_diffusion_cache/training_final/train_with_genes.csv \
  --output-dir outputs/poc_split_bundles

python scripts/leakage_and_novelty_audit.py \
  --train outputs/poc_split_bundles/random_stratified/train.txt \
  --val outputs/poc_split_bundles/random_stratified/val.txt \
  --test outputs/poc_split_bundles/random_stratified/test.txt \
  --output reports/cellline_poc_leakage_audit.tsv
```

Patch a DNA-Diffusion checkout for pretrained fine-tuning:

```bash
python scripts/patch_dnadiffusion_pretrained_finetune.py \
  --dnadiffusion-repo ~/DNA-Diffusion \
  --dry-run

python scripts/patch_dnadiffusion_pretrained_finetune.py \
  --dnadiffusion-repo ~/DNA-Diffusion
```

Outputs:

- `reports/cellline_poc_report.md`
- `reports/cellline_poc_metrics.csv`
- `reports/cellline_poc_metrics.json`
- optional figures in `reports/figures/`

The current local machine may not have the required GPU runtime. If fine-tuned generated sequences are absent, the validation report will mark the POC as incomplete and still summarize pretrained/endogenous controls.

## Stage 2: HSC/T-NK CAR Promoter System

Create a public-data manifest:

```bash
python scripts/catalog_public_hematopoietic_data.py
python scripts/curate_stage2_inputs_from_web.py --live-encode --limit-per-encode-query 25
```

Build real publishable inputs from public files:

```bash
python scripts/download_encode_stage2_public_data.py \
  --refresh-query \
  --download \
  --output-dir downloads/encode_stage2

python scripts/download_hca_stage2_expression.py \
  --output data/hsc_tnk_real/hca_expression_long.tsv

# Optional lightweight fallback if the HCA H5AD download is unavailable.
python scripts/download_hpa_stage2_expression.py \
  --output data/hsc_tnk_real/hpa_expression_long.tsv

python scripts/build_real_stage2_inputs.py \
  --output-dir data/hsc_tnk_real \
  --download-references \
  --signal-manifest downloads/encode_stage2/signal_manifest.tsv \
  --expression-long data/hsc_tnk_real/hca_expression_long.tsv \
  --tss-activity-long /path/to/tss_activity_long.tsv \
  --strict
```

The full runner does this automatically with:

```bash
STAGE2_MODE=real \
BUILD_REAL_STAGE2_INPUTS=1 \
BUILD_ENCODE_STAGE2_DATA=1 \
BUILD_HCA_EXPRESSION=1 \
BUILD_HPA_EXPRESSION=1 \
DOWNLOAD_ENCODE_FILES=1 \
DOWNLOAD_REFERENCES=1 \
bash run_all_computational.sh
```

This writes:

- `data/hsc_tnk_real/promoter_windows.tsv`
- `data/hsc_tnk_real/activity_long.tsv`
- `data/hsc_tnk_real/source_manifest.lock.tsv`
- `data/hsc_tnk_real/sample_harmonization.tsv`
- `reports/publishable/data_qc_report.md`

Merge accessibility, initiation, and expression signals:

```bash
python scripts/merge_promoter_activity_matrix.py \
  --promoters data/hsc_tnk_real/promoter_windows.tsv \
  --activity-long data/hsc_tnk_real/activity_long.tsv \
  --output data/hsc_tnk/promoter_activity_matrix.tsv
```

Select T/NK-specific promoters:

```bash
python scripts/select_tnk_specific_training_set.py \
  --activity-matrix data/hsc_tnk/promoter_activity_matrix.tsv \
  --output data/hsc_tnk/tnk_specific_training_set.tsv
```

Prepare multimodal DNA-Diffusion training data:

```bash
python scripts/train_multimodal_dnadiffusion.py \
  --training-set data/hsc_tnk/tnk_specific_training_set.tsv \
  --dnadiffusion-repo ~/DNA-Diffusion
```

Optionally patch a DNA-Diffusion checkout for true continuous profile conditioning:

```bash
python scripts/patch_dnadiffusion_multimodal_conditioning.py \
  --dnadiffusion-repo ~/DNA-Diffusion \
  --profile-dim <number_of_activity_profile_columns> \
  --dry-run

python scripts/patch_dnadiffusion_multimodal_conditioning.py \
  --dnadiffusion-repo ~/DNA-Diffusion \
  --profile-dim <number_of_activity_profile_columns>
```

Generate, rank, and design MPRA candidates:

```bash
python scripts/generate_tnk_promoters.py \
  --checkpoint /path/to/checkpoint.pt \
  --saved-data-path /path/to/hsc_tnk_multimodal_proxy.pkl \
  --output-dir outputs/hsc_tnk/generated

python scripts/filter_and_rank_candidates.py \
  --candidates outputs/hsc_tnk/generated/*.txt \
  --reference-sequences data/hsc_tnk_real/promoter_windows.tsv data/hsc_tnk/tnk_specific_training_set.tsv \
  --output outputs/hsc_tnk/ranked_candidates.tsv

python scripts/scan_candidate_motifs.py \
  --sequences outputs/hsc_tnk/ranked_candidates.tsv \
  --motifs /path/to/JASPAR_OR_HOCOMOCO.pfm \
  --output-prefix outputs/hsc_tnk/motif_validation/ranked_candidates

python scripts/annotate_promoter_ccres.py \
  --promoters data/hsc_tnk_real/promoter_windows.tsv \
  --ccre-bed /path/to/encode_screen_ccres.bed \
  --output outputs/hsc_tnk/promoter_ccre_annotations.tsv

python scripts/design_mpra_library.py \
  --ranked-candidates outputs/hsc_tnk/ranked_candidates.tsv \
  --output-prefix outputs/hsc_tnk/mpra_tnk_promoters

python scripts/validate_publishable_package.py --require-stage1
```

Analyze MPRA barcode counts after wet-lab validation:

```bash
python scripts/analyze_mpra_barcodes.py \
  --counts wetlab/barcode_counts.tsv \
  --library outputs/hsc_tnk/mpra_tnk_promoters.library.tsv \
  --output-dir outputs/hsc_tnk/mpra_analysis
```

Use MPRA results for active learning:

```bash
python scripts/design_active_learning_round.py \
  --mpra-hits outputs/hsc_tnk/mpra_analysis/sequence_level_hits.tsv \
  --candidate-features outputs/hsc_tnk/ranked_candidates.tsv \
  --output-dir outputs/hsc_tnk/active_learning_round2
```

## Important Limitation

The current upstream DNA-Diffusion code is class-label conditioned. `scripts/train_multimodal_dnadiffusion.py` preserves continuous promoter activity profiles as sidecar files and creates a profile-binned proxy TAG for immediate compatibility. `scripts/patch_dnadiffusion_multimodal_conditioning.py` applies the first-pass true-conditioning patch to a chosen checkout. The intended conditioning form is:

```text
conditioning = label_embedding(TAG) + profile_mlp(accessibility, initiation, expression)
```

in the UNet time-conditioning pathway.
