#!/bin/bash
set -euo pipefail

CHR="chr1"
WINDOWS_FILE="/path/to/Corigami_analysis/Corigami_run_dir/per_Chr_windows/chr1.txt"
OUT_ROOT="/path/to/Corigami_analysis/Corigami_run_dir"

SEQ_FILE="/path/to/corigami_data/data/hg38/dna_sequence"
CTCF_FILE="/path/to/Corigami_analysis/CTCF_REP1.mLb.clN.bigWig"
ATAC_FILE="/path/to/Corigami_analysis/24h_DMSO_merged_avg.bigWig"
MODEL_FILE="/path/to/corigami_data/model_weights/corigami_base.ckpt"
CELLTYPE="CellTypeOfInterest"

mkdir -p "$OUT_ROOT"

while read START; do
    OUTDIR="$OUT_ROOT/$CHR-$START"
    mkdir -p "$OUTDIR"

    [ -f "$OUTDIR/prediction_done.txt" ] && continue

    corigami-predict \
        --chr $CHR \
        --start $START \
        --seq $SEQ_FILE \
        --ctcf $CTCF_FILE \
        --atac $ATAC_FILE \
        --model $MODEL_FILE \
        --celltype $CELLTYPE \
        --out "$OUTDIR"

    touch "$OUTDIR/prediction_done.txt"
done < "$WINDOWS_FILE"
