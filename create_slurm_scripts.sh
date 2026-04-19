#!/bin/bash
set -euo pipefail
shopt -s nullglob

CORIGAMI_SCRIPT_DIR="/path/to/Corigami/bash_script_per_chr"
SLURM_TEMPLATE="/path/to/Corigami/Corigami_predict_slurm_submit.wilkes3"
SLURM_OUT_DIR="/path/to/Corigami/slurm_script_per_Chr"

mkdir -p "$SLURM_OUT_DIR/logs"

for SCRIPT in "$CORIGAMI_SCRIPT_DIR"/run_corigami_chr*.sh; do
    CHR=$(basename "$SCRIPT" .sh | sed 's/run_corigami_//')
    SLURM_SCRIPT="$SLURM_OUT_DIR/submit_${CHR}.slurm"

    sed \
      -e "s|^#SBATCH -J .*|#SBATCH -J Corigami_training_${CHR}|" \
      -e "s|^application=.*|application=\"$SCRIPT\"|" \
      "$SLURM_TEMPLATE" > "$SLURM_SCRIPT"

    echo "Created $SLURM_SCRIPT"
done

