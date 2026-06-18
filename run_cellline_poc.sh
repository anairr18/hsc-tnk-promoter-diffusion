#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")" && pwd)}"
DNA_DIFFUSION_REPO="${DNA_DIFFUSION_REPO:-$HOME/DNA-Diffusion}"
DNA_DIFFUSION_CACHE="${DNA_DIFFUSION_CACHE:-$HOME/dna_diffusion_cache}"
PATCH_PRETRAINED="${PATCH_PRETRAINED:-0}"
RUN_TRAINING="${RUN_TRAINING:-0}"

export DNA_DIFFUSION_REPO DNA_DIFFUSION_CACHE

python "$PROJECT_ROOT/scripts/download_project_assets.py" --output-dir "$PROJECT_ROOT/downloads"
python "$PROJECT_ROOT/scripts/check_pretrained_initialization.py" \
  --dnadiffusion-repo "$DNA_DIFFUSION_REPO" \
  --output "$PROJECT_ROOT/reports/pretrained_initialization_check.json"

if [[ "$PATCH_PRETRAINED" == "1" ]]; then
  python "$PROJECT_ROOT/scripts/patch_dnadiffusion_pretrained_finetune.py" \
    --dnadiffusion-repo "$DNA_DIFFUSION_REPO"
fi

if [[ "$RUN_TRAINING" == "1" ]]; then
  python "$PROJECT_ROOT/pipeline/cell1_dataset_and_training.py"
  python "$PROJECT_ROOT/pipeline/cell2_generate_and_backup.py"
fi

python "$PROJECT_ROOT/scripts/validate_cellline_poc.py" --project-root "$PROJECT_ROOT"

SPLIT_DIR="$DNA_DIFFUSION_CACHE/training_final"
GENERATED_DIR="$SPLIT_DIR/generated"
if [[ -f "$SPLIT_DIR/train.txt" && -f "$SPLIT_DIR/test.txt" ]]; then
  mapfile -t GENERATED_FILES < <(find "$GENERATED_DIR" -maxdepth 1 -type f 2>/dev/null || true)
  python "$PROJECT_ROOT/scripts/leakage_and_novelty_audit.py" \
    --train "$SPLIT_DIR/train.txt" \
    --val "$SPLIT_DIR/val.txt" \
    --test "$SPLIT_DIR/test.txt" \
    --generated "${GENERATED_FILES[@]}" \
    --output "$PROJECT_ROOT/reports/cellline_poc_leakage_audit.tsv"
fi

echo "Cell-line POC workflow complete. Read reports/cellline_poc_report.md"
