#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")" && pwd)}"
BASE_DIR="${BASE_DIR:-$(dirname "$PROJECT_ROOT")}"
DNA_DIFFUSION_REPO="${DNA_DIFFUSION_REPO:-$BASE_DIR/DNA-Diffusion}"
DNA_DIFFUSION_CACHE="${DNA_DIFFUSION_CACHE:-$BASE_DIR/dna_diffusion_cache}"

RUN_STAGE1="${RUN_STAGE1:-1}"
RUN_STAGE2="${RUN_STAGE2:-1}"
STAGE2_MODE="${STAGE2_MODE:-demo}"   # demo or real
BUILD_REAL_STAGE2_INPUTS="${BUILD_REAL_STAGE2_INPUTS:-0}"
STRICT_REAL_INPUTS="${STRICT_REAL_INPUTS:-1}"
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
  if [[ "$BUILD_REAL_STAGE2_INPUTS" == "1" ]]; then
    REAL_INPUT_DIR="${REAL_INPUT_DIR:-$PROJECT_ROOT/data/hsc_tnk_real}"
    REAL_BUILD_ARGS=(
      "$PROJECT_ROOT/scripts/build_real_stage2_inputs.py"
      --output-dir "$REAL_INPUT_DIR"
    )
    if [[ "${DOWNLOAD_REFERENCES:-0}" == "1" ]]; then
      REAL_BUILD_ARGS+=(--download-references)
    fi
    if [[ -n "${HG38_FASTA:-}" ]]; then
      REAL_BUILD_ARGS+=(--hg38-fasta "$HG38_FASTA")
    fi
    if [[ -n "${GENCODE_GTF:-}" ]]; then
      REAL_BUILD_ARGS+=(--gencode-gtf "$GENCODE_GTF")
    fi
    if [[ -n "${FANTOM_TSS_BED:-}" ]]; then
      REAL_BUILD_ARGS+=(--fantom-tss-bed "$FANTOM_TSS_BED")
    fi
    if [[ -n "${SIGNAL_MANIFEST:-}" ]]; then
      REAL_BUILD_ARGS+=(--signal-manifest "$SIGNAL_MANIFEST")
    fi
    if [[ -n "${EXPRESSION_LONG_TSV:-}" ]]; then
      REAL_BUILD_ARGS+=(--expression-long "$EXPRESSION_LONG_TSV")
    fi
    if [[ -n "${TSS_ACTIVITY_LONG_TSV:-}" ]]; then
      REAL_BUILD_ARGS+=(--tss-activity-long "$TSS_ACTIVITY_LONG_TSV")
    fi
    if [[ -n "${MAX_PROMOTERS:-}" ]]; then
      REAL_BUILD_ARGS+=(--max-promoters "$MAX_PROMOTERS")
    fi
    if [[ "$STRICT_REAL_INPUTS" == "1" ]]; then
      REAL_BUILD_ARGS+=(--strict)
    fi
    python "${REAL_BUILD_ARGS[@]}"
    PROMOTER_WINDOWS="$REAL_INPUT_DIR/promoter_windows.tsv"
    ACTIVITY_LONG_TSV="$REAL_INPUT_DIR/activity_long.tsv"
  elif [[ -n "${PROMOTER_WINDOWS:-}" && -n "${ACTIVITY_LONG_TSV:-}" ]]; then
    PROMOTER_WINDOWS="$PROMOTER_WINDOWS"
    ACTIVITY_LONG_TSV="$ACTIVITY_LONG_TSV"
  else
    echo "Real Stage 2 requires either BUILD_REAL_STAGE2_INPUTS=1 or both PROMOTER_WINDOWS and ACTIVITY_LONG_TSV." >&2
    exit 1
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
  --reference-sequences "$PROMOTER_WINDOWS" "$STAGE2_DIR/tnk_specific_training_set.tsv" \
  --output "$STAGE2_DIR/ranked_candidates.tsv" \
  --top-n 600

if [[ -n "${MOTIF_FILES:-}" ]]; then
  # shellcheck disable=SC2206
  MOTIF_ARRAY=($MOTIF_FILES)
  python "$PROJECT_ROOT/scripts/scan_candidate_motifs.py" \
    --sequences "$STAGE2_DIR/ranked_candidates.tsv" \
    --motifs "${MOTIF_ARRAY[@]}" \
    --output-prefix "$STAGE2_DIR/motif_validation/ranked_candidates"
fi

if [[ -n "${CCRE_BED:-}" ]]; then
  python "$PROJECT_ROOT/scripts/annotate_promoter_ccres.py" \
    --promoters "$PROMOTER_WINDOWS" \
    --ccre-bed "$CCRE_BED" \
    --output "$STAGE2_DIR/promoter_ccre_annotations.tsv"
fi

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

if [[ "$STAGE2_MODE" == "real" ]]; then
  python "$PROJECT_ROOT/scripts/validate_publishable_package.py" \
    --require-stage1 \
    --promoter-windows "$PROMOTER_WINDOWS" \
    --activity-long "$ACTIVITY_LONG_TSV" \
    --ranked-candidates "$STAGE2_DIR/ranked_candidates.tsv" \
    --mpra-library "$STAGE2_DIR/mpra_tnk_promoters.library.tsv"
fi

echo "DONE. See reports/all_computational_summary.md"
