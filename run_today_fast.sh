#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")" && pwd)}"
DNA_DIFFUSION_REPO="${DNA_DIFFUSION_REPO:-$HOME/DNA-Diffusion}"
DNA_DIFFUSION_CACHE="${DNA_DIFFUSION_CACHE:-$HOME/dna_diffusion_cache}"
ALLOW_CPU="${ALLOW_CPU:-0}"
PATCH_PRETRAINED="${PATCH_PRETRAINED:-0}"
DNA_DIFFUSION_EPOCHS="${DNA_DIFFUSION_EPOCHS:-5000}"
DNA_DIFFUSION_MIN_EPOCHS="${DNA_DIFFUSION_MIN_EPOCHS:-2000}"

export DNA_DIFFUSION_REPO DNA_DIFFUSION_CACHE DNA_DIFFUSION_EPOCHS DNA_DIFFUSION_MIN_EPOCHS

echo "== Checking GPU =="
if ! python - <<'PY'
import torch
print("cuda", torch.cuda.is_available())
print("device_count", torch.cuda.device_count())
raise SystemExit(0 if torch.cuda.is_available() else 2)
PY
then
  if [[ "$ALLOW_CPU" != "1" ]]; then
    echo "No CUDA GPU detected. Re-run on a GPU machine, or set ALLOW_CPU=1 for a slow CPU smoke run." >&2
    exit 2
  fi
fi

echo "== Installing lightweight requirements =="
python -m pip install -r "$PROJECT_ROOT/requirements.txt"

echo "== Ensuring uv and DNA-Diffusion repo =="
if ! command -v uv >/dev/null 2>&1; then
  python -m pip install -q uv
fi
if [[ ! -d "$DNA_DIFFUSION_REPO/.git" ]]; then
  rm -rf "$DNA_DIFFUSION_REPO"
  git clone https://github.com/pinellolab/DNA-Diffusion.git "$DNA_DIFFUSION_REPO"
fi

echo "== Downloading metadata =="
python "$PROJECT_ROOT/scripts/download_project_assets.py" --output-dir "$PROJECT_ROOT/downloads"

echo "== Checking pretrained fine-tune support =="
python "$PROJECT_ROOT/scripts/check_pretrained_initialization.py" \
  --dnadiffusion-repo "$DNA_DIFFUSION_REPO" \
  --output "$PROJECT_ROOT/reports/pretrained_initialization_check.json"

if [[ "$PATCH_PRETRAINED" == "1" ]]; then
  echo "== Patching DNA-Diffusion train.py for pretrained fine-tuning =="
  python "$PROJECT_ROOT/scripts/patch_dnadiffusion_pretrained_finetune.py" \
    --dnadiffusion-repo "$DNA_DIFFUSION_REPO"
fi

echo "== Running cell-line POC dataset + training =="
python "$PROJECT_ROOT/pipeline/cell1_dataset_and_training.py"

echo "== Generating sequences =="
python "$PROJECT_ROOT/pipeline/cell2_generate_and_backup.py"

echo "== Validating POC =="
python "$PROJECT_ROOT/scripts/validate_cellline_poc.py" --project-root "$PROJECT_ROOT"

echo "DONE. Open reports/cellline_poc_report.md"
