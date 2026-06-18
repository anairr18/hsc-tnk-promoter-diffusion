#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "$0")" && pwd)}"
DNA_DIFFUSION_REPO="${DNA_DIFFUSION_REPO:-$HOME/DNA-Diffusion}"
DATA_DIR="$PROJECT_ROOT/data/hsc_tnk"
mkdir -p "$DATA_DIR"

python "$PROJECT_ROOT/scripts/catalog_public_hematopoietic_data.py"

if [[ -n "${HG38_FASTA:-}" && (-n "${GENCODE_GTF:-}" || -n "${FANTOM_TSS_BED:-}") ]]; then
  CMD=(python "$PROJECT_ROOT/scripts/build_promoter_windows.py" --fasta "$HG38_FASTA" --output "$DATA_DIR/promoter_windows.tsv")
  [[ -n "${GENCODE_GTF:-}" ]] && CMD+=(--gencode-gtf "$GENCODE_GTF")
  [[ -n "${FANTOM_TSS_BED:-}" ]] && CMD+=(--fantom-tss-bed "$FANTOM_TSS_BED")
  "${CMD[@]}"
fi

if [[ -n "${ACTIVITY_LONG_TSV:-}" ]]; then
  python "$PROJECT_ROOT/scripts/merge_promoter_activity_matrix.py" \
    --promoters "$DATA_DIR/promoter_windows.tsv" \
    --activity-long "$ACTIVITY_LONG_TSV" \
    --output "$DATA_DIR/promoter_activity_matrix.tsv"

  python "$PROJECT_ROOT/scripts/select_tnk_specific_training_set.py" \
    --activity-matrix "$DATA_DIR/promoter_activity_matrix.tsv" \
    --output "$DATA_DIR/tnk_specific_training_set.tsv"

  python "$PROJECT_ROOT/scripts/train_multimodal_dnadiffusion.py" \
    --training-set "$DATA_DIR/tnk_specific_training_set.tsv" \
    --dnadiffusion-repo "$DNA_DIFFUSION_REPO"
fi

echo "HSC/T-NK backbone prepared under $DATA_DIR"
