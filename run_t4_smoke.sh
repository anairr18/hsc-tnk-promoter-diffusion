#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")" && pwd)}"
BASE_DIR="${BASE_DIR:-$(dirname "$PROJECT_ROOT")}"
DNA_DIFFUSION_REPO="${DNA_DIFFUSION_REPO:-$BASE_DIR/DNA-Diffusion}"
PATCH_PRETRAINED="${PATCH_PRETRAINED:-1}"
SMOKE_EPOCHS="${SMOKE_EPOCHS:-3}"
SMOKE_MIN_EPOCHS="${SMOKE_MIN_EPOCHS:-1}"

python -m pip install -q -r "$PROJECT_ROOT/requirements.txt"

ARGS=(
  python "$PROJECT_ROOT/scripts/colab_t4_smoke.py"
  --project-root "$PROJECT_ROOT"
  --dnadiffusion-repo "$DNA_DIFFUSION_REPO"
  --epochs "$SMOKE_EPOCHS"
  --min-epochs "$SMOKE_MIN_EPOCHS"
)

if [[ "$PATCH_PRETRAINED" == "1" ]]; then
  ARGS+=(--patch-pretrained)
fi

"${ARGS[@]}"
