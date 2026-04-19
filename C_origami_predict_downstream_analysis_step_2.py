#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 16:46:30 2026

@author: d
"""

#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pyranges as pr
import os

import matplotlib.pyplot as plt

# =========================
# INPUTS
# =========================
outdir ="/path/to/Corigami/Corigami_results_summary"
figdir = "/path/to/Corigami/Corigami_results_summary/figures"

loops_file = "Corigami_predict_downstream_analysis_origami_loops.bedpe"
chip_file = "sorted_CellTypeOfInterest_GeneOfInterest_narrowpeak_consensus.bed"

out_prefix = "Corigami_GeneOfInterest_results"

# =========================
# LOAD RAW BEDPE (STRICT)
# =========================
loops_df = pd.read_csv(
    loops_file,
    sep="\t",
    header=None,
    names=[
        "chrom1","start1","end1",
        "chrom2","start2","end2",
        "score","distance"
    ]
)

print("Raw loops:", loops_df.shape)

# =========================
# CLEAN TYPES
# =========================
loops_df["chrom1"] = loops_df["chrom1"].astype(str)
loops_df["chrom2"] = loops_df["chrom2"].astype(str)

for col in ["start1","end1","start2","end2","score","distance"]:
    loops_df[col] = pd.to_numeric(loops_df[col], errors="coerce")

loops_df = loops_df.dropna(subset=["start1","start2","end1","end2","score"])

loops_df[["start1","end1","start2","end2"]] = loops_df[
    ["start1","end1","start2","end2"]
].astype(int)

print("After cleaning:", loops_df.shape)

if len(loops_df) == 0:
    raise ValueError(" No loops after cleaning")

# =========================
# LOAD GeneOfInterest PEAKS
# =========================
GeneOfInterest = pr.read_bed(chip_file)

# =========================
# BUILD PyRanges (NO loop_id!)
# =========================
anchor1 = pr.PyRanges(pd.DataFrame({
    "Chromosome": loops_df["chrom1"],
    "Start": loops_df["start1"],
    "End": loops_df["end1"]
}))

anchor2 = pr.PyRanges(pd.DataFrame({
    "Chromosome": loops_df["chrom2"],
    "Start": loops_df["start2"],
    "End": loops_df["end2"]
}))

# =========================
# INTERSECTIONS
# =========================
a1 = anchor1.join(GeneOfInterest)
a2 = anchor2.join(GeneOfInterest)

# match by row index (safe + stable)
a1_hits = set(a1.df.index)
a2_hits = set(a2.df.index)

# =========================
# ANNOTATE IN MEMORY ONLY
# =========================
loops_df["anchor1_GeneOfInterest"] = loops_df.index.isin(a1_hits)
loops_df["anchor2_GeneOfInterest"] = loops_df.index.isin(a2_hits)

loops_df["any_GeneOfInterest"] = loops_df["anchor1_GeneOfInterest"] | loops_df["anchor2_GeneOfInterest"]
loops_df["both_GeneOfInterest"] = loops_df["anchor1_GeneOfInterest"] & loops_df["anchor2_GeneOfInterest"]

# =========================
# SUMMARY STATS
# =========================
print("\n========== GeneOfInterest SUMMARY ==========")
print("Total loops:", len(loops_df))
print("GeneOfInterest loops:", loops_df["any_GeneOfInterest"].sum())
print("Fraction GeneOfInterest:", loops_df["any_GeneOfInterest"].mean())
print("GeneOfInterest-GeneOfInterest loops:", loops_df["both_GeneOfInterest"].sum())

print("\nMean loop score:")
print(loops_df.groupby("any_GeneOfInterest")["score"].mean())



# ========== GeneOfInterest SUMMARY ==========
# Total loops: 12182
# GeneOfInterest loops: 100
# Fraction GeneOfInterest: 0.008208832703989493
# GeneOfInterest-GeneOfInterest loops: 79

# Mean loop score:
# any_GeneOfInterest
# False    7.975422
# True     8.680740
# Name: score, dtype: float64

# =========================
# EXPORT CLEAN BEDPE (IMPORTANT)
# =========================
bedpe_out = f"{out_prefix}_clean_loops.bedpe"

loops_df[[
    "chrom1","start1","end1",
    "chrom2","start2","end2",
    "score","distance"
]].to_csv(
    bedpe_out,
    sep="\t",
    header=False,
    index=False
)

print("\nSaved clean BEDPE:", bedpe_out)

# =========================
# EXPORT ANNOTATIONS (SEPARATE FILE)
# =========================
annot_out = f"{out_prefix}_GeneOfInterest_annotations.tsv"

loops_df.to_csv(
    annot_out,
    sep="\t",
    index=False
)

print("Saved annotation table:", annot_out)


# =========================
# FIGURE 1: GeneOfInterest FRACTION
# =========================
plt.figure(figsize=(4,4))
plt.bar(["GeneOfInterest loops", "Other loops"],
        [loops_df["any_GeneOfInterest"].mean(),
         1 - loops_df["any_GeneOfInterest"].mean()])
plt.ylabel("Fraction of loops")
plt.title("GeneOfInterest association with loops")
plt.tight_layout()
plt.savefig(os.path.join(figdir, "GeneOfInterest_loop_fraction.png"), dpi=500)
plt.close()


# =========================
# FIGURE 2: LOOP DISTANCE DISTRIBUTION
# =========================
plt.figure(figsize=(5,4))
plt.hist(loops_df["distance"], bins=50)
plt.xlabel("Loop distance")
plt.ylabel("Count")
plt.title("Predicted loop distances")
plt.tight_layout()
plt.savefig(os.path.join(figdir, "loop_distance_distribution.png"), dpi=500)
plt.close()


# =========================
# FIGURE 3: TOP LOOP REGIONS (SIMPLE SNAPSHOT LIST)
# =========================
top_loops = loops_df.sort_values("score", ascending=False).head(5)

top_loops.to_csv(os.path.join(outdir, "top_loops.tsv"), sep="\t", index=False)


# =========================
# EXPORT FINAL TABLES
# =========================
loops_df.to_csv(os.path.join(outdir, "Corigami_loops_annotated.tsv"),
                sep="\t", index=False)

print("\n================ DONE ================\n")
print("Output folder:", outdir)
print("Key files:")
print("- GeneOfInterest_summary.txt")
print("- Corigami_loops_annotated.tsv")
print("- figures/")
print("=====================================\n")


# =========================
# IGV EXPORT
# =========================

track_dir = os.path.join(outdir, "IGV_tracks")
os.makedirs(track_dir, exist_ok=True)

# -------------------------
# 1. Loop BEDPE 
# -------------------------
loops_df[[
    "chrom1","start1","end1",
    "chrom2","start2","end2","score"
]].to_csv(
    os.path.join(track_dir, "Corigami_loops.bedpe"),
    sep="\t", header=False, index=False
)

# -------------------------
# 2. GeneOfInterest peaks (clean copy for IGV)
# -------------------------
GeneOfInterest.df[["Chromosome", "Start", "End"]].to_csv(
    os.path.join(track_dir, "GeneOfInterest_peaks.bed"),
    sep="\t",
    header=False,
    index=False
)

# -------------------------
# 3. Simple BEDGRAPH-style loop density track
# -------------------------
# (counts loop anchors per region bin)

bin_size = 50000  # coarse genome bins for visualization

all_positions = []

for _, row in loops_df.iterrows():
    all_positions.append((row["chrom1"], row["start1"]))
    all_positions.append((row["chrom2"], row["start2"]))

pos_df = pd.DataFrame(all_positions, columns=["chrom", "pos"])

pos_df["bin"] = (pos_df["pos"] // bin_size) * bin_size

track = pos_df.groupby(["chrom", "bin"]).size().reset_index(name="value")

track.to_csv(
    os.path.join(track_dir, "loop_anchor_density.bedgraph"),
    sep="\t", header=False, index=False
)