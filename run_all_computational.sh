#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")" && pwd)}"
BASE_DIR="${BASE_DIR:-$(dirname "$PROJECT_ROOT")}"
DNA_DIFFUSION_REPO="${DNA_DIFFUSION_REPO:-$BASE_DIR/DNA-Diffusion}"
DNA_DIFFUSION_CACHE="${DNA_DIFFUSION_CACHE:-$BASE_DIR/dna_diffusion_cache}"

RUN_STAGE1="${RUN_STAGE1:-1}"
RUN_STAGE2="${RUN_STAGE2:-1}"
STAGE2_MODE="${STAGE2_MODE:-demo}"   # demo or real
PATCH_PRETRAINED="${PATCH_PRETRAINED:-1}"

DNA_DIFFUSION_EPOCHS="${DNA_DIFFUSION_EPOCHS:-3000}"
DNA_DIFFUSION_MIN_EPOCHS="${DNA_DIFFUSION_MIN_EPOCHS:-1000}"
STAGE2_EPOCHS="${STAGE2_EPOCHS:-500}"
STAGE2_MIN_EPOCHS="${STAGE2_MIN_EPOCHS:-100}"
STAGE2_BATCH_SIZE="${STAGE2_BATCH_SIZE:-64}"
STAGE2_SAMPLES="${STAGE2_SAMPLES:-1000}"

export DNA_DIFFUSION_REPO DNA_DIFFUSION_CACHE PATCH_PRETRAINED DNA_DIFFUSION_EPOCHS DNA_DIFFUSION_MIN_EPOCHS

echo "== Installing requirements =="
python -m pip install -q -r "$PROJECT_ROOT/requirements.txt"

echo "== Stage 0: manifests/assets =="
python "$PROJECT_ROOT/scripts/download_project_assets.py" --output-dir "$PROJECT_ROOT/downloads"
python "$PROJECT_ROOT/scripts/catalog_public_hematopoietic_data.py"
python "$PROJECT_ROOT/scripts/curate_stage2_inputs_from_web.py" --live-encode --limit-per-encode-query 25

if [[ "$RUN_STAGE1" == "1" ]]; then
  echo "== Stage 1: cell-line POC =="
  bash "$PROJECT_ROOT/run_today_fast.sh"
fi

if [[ "$RUN_STAGE2" != "1" ]]; then
  echo "Stage 2 skipped."
  exit 0
fi

echo "== Stage 2: HSC/T-NK CAR promoter computational design =="
STAGE2_DIR="$PROJECT_ROOT/outputs/hsc_tnk"
INPUT_DIR="$PROJECT_ROOT/data/hsc_tnk"
mkdir -p "$STAGE2_DIR" "$INPUT_DIR"

if [[ "$STAGE2_MODE" == "demo" ]]; then
  echo "Using demo Stage 2 inputs. This validates the full pipeline but is not biological evidence."
  python "$PROJECT_ROOT/scripts/build_demo_hsc_tnk_inputs.py" --output-dir "$INPUT_DIR/demo_inputs"
  PROMOTER_WINDOWS="$INPUT_DIR/demo_inputs/promoter_windows.tsv"
  ACTIVITY_LONG_TSV="$INPUT_DIR/demo_inputs/activity_long.tsv"
else
  : "${ACTIVITY_LONG_TSV:?Set ACTIVITY_LONG_TSV for real Stage 2 mode}"
  if [[ -n "${PROMOTER_WINDOWS:-}" ]]; then
    PROMOTER_WINDOWS="$PROMOTER_WINDOWS"
  else
    : "${HG38_FASTA:?Set HG38_FASTA when PROMOTER_WINDOWS is not supplied}"
    : "${GENCODE_GTF:?Set GENCODE_GTF or provide PROMOTER_WINDOWS}"
    PROMOTER_WINDOWS="$INPUT_DIR/promoter_windows.tsv"
    python "$PROJECT_ROOT/scripts/build_promoter_windows.py" \
      --fasta "$HG38_FASTA" \
      --gencode-gtf "$GENCODE_GTF" \
      ${FANTOM_TSS_BED:+--fantom-tss-bed "$FANTOM_TSS_BED"} \
      --output "$PROMOTER_WINDOWS"
  fi
fi

python "$PROJECT_ROOT/scripts/merge_promoter_activity_matrix.py" \
  --promoters "$PROMOTER_WINDOWS" \
  --activity-long "$ACTIVITY_LONG_TSV" \
  --output "$STAGE2_DIR/promoter_activity_matrix.tsv"

python "$PROJECT_ROOT/scripts/select_tnk_specific_training_set.py" \
  --activity-matrix "$STAGE2_DIR/promoter_activity_matrix.tsv" \
  --output "$STAGE2_DIR/tnk_specific_training_set.tsv"

python "$PROJECT_ROOT/scripts/train_multimodal_dnadiffusion.py" \
  --training-set "$STAGE2_DIR/tnk_specific_training_set.tsv" \
  --dnadiffusion-repo "$DNA_DIFFUSION_REPO" \
  --output-dir "$STAGE2_DIR/multimodal_training" \
  --run \
  --epochs "$STAGE2_EPOCHS" \
  --min-epochs "$STAGE2_MIN_EPOCHS" \
  --batch-size "$STAGE2_BATCH_SIZE"

LATEST_CKPT="$(ls -t "$DNA_DIFFUSION_REPO"/checkpoints/*.pt | head -n 1)"
SAVED_DATA="$STAGE2_DIR/multimodal_training/hsc_tnk_multimodal_proxy.pkl"

python "$PROJECT_ROOT/scripts/generate_tnk_promoters.py" \
  --dnadiffusion-repo "$DNA_DIFFUSION_REPO" \
  --checkpoint "$LATEST_CKPT" \
  --saved-data-path "$SAVED_DATA" \
  --cell-types "TNK_HIGH_PROFILE_Q3" \
  --number-of-samples "$STAGE2_SAMPLES" \
  --output-dir "$STAGE2_DIR/generated" \
  --run

python "$PROJECT_ROOT/scripts/filter_and_rank_candidates.py" \
  --candidates "$STAGE2_DIR"/generated/*.txt \
  --output "$STAGE2_DIR/ranked_candidates.tsv" \
  --top-n 600

python "$PROJECT_ROOT/scripts/create_mpra_control_set.py" \
  --activity-matrix "$STAGE2_DIR/promoter_activity_matrix.tsv" \
  --output-dir "$STAGE2_DIR/mpra_controls" \
  --n 96

python "$PROJECT_ROOT/scripts/design_mpra_library.py" \
  --ranked-candidates "$STAGE2_DIR/ranked_candidates.tsv" \
  --control-files "$STAGE2_DIR"/mpra_controls/*.txt \
  --output-prefix "$STAGE2_DIR/mpra_tnk_promoters" \
  --barcodes-per-sequence 12

cat > "$PROJECT_ROOT/reports/all_computational_summary.md" <<EOF
# All Computational Workflow Summary

- Stage 1 POC report: \`reports/cellline_poc_report.md\`
- Stage 2 mode: \`$STAGE2_MODE\`
- Stage 2 activity matrix: \`outputs/hsc_tnk/promoter_activity_matrix.tsv\`
- Stage 2 selected training set: \`outputs/hsc_tnk/tnk_specific_training_set.tsv\`
- Stage 2 generated candidates: \`outputs/hsc_tnk/generated/\`
- Stage 2 ranked candidates: \`outputs/hsc_tnk/ranked_candidates.tsv\`
- MPRA library design: \`outputs/hsc_tnk/mpra_tnk_promoters.library.tsv\`

Note: \`demo\` Stage 2 mode validates software execution but does not constitute biological evidence. Use \`STAGE2_MODE=real\` with real hematopoietic activity inputs for real results.
EOF

echo "DONE. See reports/all_computational_summary.md"
