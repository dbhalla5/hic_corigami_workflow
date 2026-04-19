#!/bin/bash

# ===== HPC PATHS (LITERAL STRINGS) =====
HPC_BASE="/path/to/Corigami_analysis/Corigami_run_dir"
HPC_WINDOWS_DIR="$HPC_BASE/per_Chr_windows"
HPC_OUT_DIR="$HPC_BASE"

# ===== LOCAL PATH =====
BASH_SCRIPTS_DIR="/path/to/Corigami_analysis/for_HPC/bash_script_per_chr"
mkdir -p "$BASH_SCRIPTS_DIR"

# Corigami inputs (HPC)
SEQ_FILE="/path/to/corigami_data/data/hg38/dna_sequence"
CTCF_FILE="/path/to/Corigami_analysis/CTCF_REP1.mLb.clN.bigWig"
ATAC_FILE="/path/to/Corigami_analysis/24h_DMSO_merged_avg.bigWig"
MODEL_FILE="/path/to/corigami_data/model_weights/corigami_base.ckpt"
CELLTYPE="CellTypeOfInterest"

# ===== Chromosome list (explicit, no FS access) =====
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

for CHR in "${CHROMS[@]}"; do
    SCRIPT_FILE="$BASH_SCRIPTS_DIR/run_corigami_${CHR}.sh"

    cat > "$SCRIPT_FILE" <<EOF
#!/bin/bash
set -euo pipefail

CHR="$CHR"
WINDOWS_FILE="$HPC_WINDOWS_DIR/${CHR}.txt"
OUT_ROOT="$HPC_OUT_DIR"

SEQ_FILE="$SEQ_FILE"
CTCF_FILE="$CTCF_FILE"
ATAC_FILE="$ATAC_FILE"
MODEL_FILE="$MODEL_FILE"
CELLTYPE="$CELLTYPE"

mkdir -p "\$OUT_ROOT"

while read _ START; do
    OUTDIR="\$OUT_ROOT/\$CHR-\$START"
    mkdir -p "\$OUTDIR"

    [ -f "\$OUTDIR/prediction_done.txt" ] && continue

    corigami-predict \\
        --chr \$CHR \\
        --start \$START \\
        --seq \$SEQ_FILE \\
        --ctcf \$CTCF_FILE \\
        --atac \$ATAC_FILE \\
        --model \$MODEL_FILE \\
        --celltype \$CELLTYPE \\
        --out "\$OUTDIR"

    touch "\$OUTDIR/prediction_done.txt"
done < "\$WINDOWS_FILE"
EOF

    chmod +x "$SCRIPT_FILE"
    echo "Created: $SCRIPT_FILE"
done
