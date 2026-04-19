#!/bin/bash

# -----------------------------
# User-editable paths
# -----------------------------
# Path to the master windows file (all chromosomes)
WINDOWS_FILE="/path/to/for_HPC/corigami_windows_canonical.txt"

# Root folder where per-chromosome folders will be written
OUTPUT_ROOT="/path/to/for_HPC/per_Chr_windows"



# -----------------------------
# Prepare output directories
# -----------------------------
mkdir -p "$OUTPUT_ROOT"


# -----------------------------
# Split master file into per-chromosome files
# -----------------------------
# Master file expected format: CHR<TAB>START
# e.g., chr1   4943993
# ------------------------------------------------
echo "Splitting $WINDOWS_FILE into per-chromosome files..."

# Get all chromosomes in the master file
CHRS=$(awk '{print $1}' "$WINDOWS_FILE" | sort | uniq)

for CHR in $CHRS; do
    # File to hold start positions for this chromosome
    CHR_FILE="$OUTPUT_ROOT/${CHR}.txt"

    # Extract all start positions for this chromosome
    awk -v chr="$CHR" '$1==chr {print $2}' "$WINDOWS_FILE" > "$CHR_FILE"

    echo "  $CHR -> $CHR_FILE ($(wc -l < "$CHR_FILE") windows)"
done

echo "Per-chromosome window files ready in $OUTPUT_ROOT"

