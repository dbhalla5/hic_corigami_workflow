#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 14:08:58 2026

@author: d
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import maximum_filter

import pyranges as pr
from scipy.stats import mannwhitneyu
from scipy.stats import zscore

root = "/path/to/Corigami/prediction_results"   # folder containing chr1, chr2, ...

chip_file = "sorted_CellTypeOfInterest_GeneOfInterest_narrowpeak_consensus.bed"
bed_file = "Corigami_predict_downstream_analysis_CellTypeOfInterest_TAD_boundaries.bed"
bedpe_file = "Corigami_predict_downstream_analysis_origami_loops.bedpe"


records = []

for chr_dir in sorted(os.listdir(root)):
    if not chr_dir.startswith("chr"):
        continue

    chr_path = os.path.join(root, chr_dir)
    if not os.path.isdir(chr_path):
        continue

    for window_dir in os.listdir(chr_path):
        if "-" not in window_dir:
            continue

        chr_name, start_str = window_dir.split("-")
        start = int(start_str)

        npy_path = os.path.join(
            chr_path, window_dir, "CellTypeOfInterest", "prediction", "npy",
            f"{chr_name}_{start}.npy"
        )

        if not os.path.exists(npy_path):
            continue

        arr = np.load(npy_path)

        records.append({
            "chromosome": chr_name,
            "start": start,
            "path": npy_path,
            "mean": float(arr.mean()),
            "std": float(arr.std()),
            "shape": arr.shape
        })

df = pd.DataFrame(records).sort_values(["chromosome", "start"])
df


##  Visualize the predicted contact maps

def plot_window(path, vmin=-0.1, vmax=0.1):
    arr = np.load(path)
    plt.figure(figsize=(5,5))
    plt.imshow(arr, cmap="RdBu_r", vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.title(os.path.basename(path))
    plt.show()
    
    
# 1. Insulation score (simple implementation)


def compute_insulation(matrix, window_bins=5):
    """
    Simple insulation score along the diagonal.
    window_bins: half-size of the square around each diagonal bin.
    """
    n = matrix.shape[0]
    ins = np.full(n, np.nan)

    for i in range(window_bins, n - window_bins):
        # square around the diagonal position (i, i)
        sub = matrix[i - window_bins:i + window_bins + 1,
                     i - window_bins:i + window_bins + 1]
        ins[i] = np.nanmean(sub)

    return ins


# 2. loop detection 

def distance_normalize(matrix):
    """
    Normalize contact matrix by genomic distance (diagonal-wise).
    Removes distance decay effect.
    """
    mat = matrix.copy()
    n = mat.shape[0]

    for d in range(n):
        diag = np.diag(mat, k=d)
        if len(diag) == 0:
            continue

        mean_val = np.nanmean(diag)
        if mean_val == 0 or np.isnan(mean_val):
            continue

        i = np.arange(n - d)
        j = np.arange(d, n)

        mat[i, j] /= mean_val
        mat[j, i] /= mean_val  # symmetry

    return mat

def detect_loops(matrix, neighborhood=3, z_thresh=3.0, min_dist=10):
    """
    Improved loop detection:
    - distance-normalized input
    - stricter threshold
    - removes short-range artifacts
    """
    arr = matrix.copy()

    # z-score normalize
    m = np.nanmean(arr)
    s = np.nanstd(arr)
    if s == 0 or np.isnan(s):
        return []

    z = (arr - m) / s

    # local maxima
    max_f = maximum_filter(z, size=neighborhood, mode='constant', cval=-np.inf)
    peaks = (z == max_f) & (z >= z_thresh)

    coords = np.argwhere(peaks)

    # filter by genomic distance (remove diagonal noise)
    loops = []


    for i, j in coords:
        
        if abs(i - j) < min_dist:
            continue
        if i >= j:   # ← ADD THIS LINE
            continue
        loops.append((int(i), int(j), float(z[i, j])))

    return loops


#######################


#  3. Apply to all windows in the dataframe

def analyze_window(path):
    arr = np.load(path)

    ins = compute_insulation(arr, window_bins=5)
    # simple boundary strength: inverse of insulation (lower insulation → stronger boundary)
    # can also look at local minima explicitly
    boundary_strength = float(np.nanstd(ins))  # crude global measure


  
    arr_norm = distance_normalize(arr)
    loops = detect_loops(arr_norm, neighborhood=3, z_thresh=3.0)
    
    n_loops = len(loops)

    return ins, boundary_strength, n_loops, loops

results = []

for idx, row in df.iterrows():
    arr = np.load(row["path"])
    ins = compute_insulation(arr, window_bins=5)
    boundary_strength = float(np.nanstd(ins))
    
    
    arr_norm = distance_normalize(arr)
    loops = detect_loops(arr_norm, neighborhood=3, z_thresh=3.0)   
     
    n_loops = len(loops)

    results.append({
        "chromosome": row["chromosome"],
        "start": row["start"],
        "path": row["path"],
        "mean": row["mean"],
        "std": row["std"],
        "boundary_strength": boundary_strength,
        "n_loops": n_loops,
    })

df_struct = pd.DataFrame(results).sort_values(["chromosome", "start"])
df_struct.head()


##############################################################################

# Build genome‑wide insulation tracks

all_insulation = []

for idx, row in df.iterrows():
    arr = np.load(row["path"])
    ins = compute_insulation(arr, window_bins=5)
    # convert bin index → genomic coordinate
    bin_size = 2048  # C.Origami default
    coords = row["start"] + np.arange(len(ins)) * bin_size
    for c, v in zip(coords, ins):
        all_insulation.append([row["chromosome"], c, v])

ins_df = pd.DataFrame(all_insulation, columns=["chrom", "pos", "insulation"])


# Call TAD boundaries

ins_df["d_ins"] = ins_df.groupby("chrom")["insulation"].diff()
boundaries = ins_df[(ins_df["d_ins"] < 0) & (ins_df["d_ins"].shift(-1) > 0)]


# Export boundaries as BED

boundaries[["chrom", "pos", "pos"]].to_csv( bed_file ,
                                           sep="\t", header=False, index=False)

############################################################################


# Map structures to genomic coordinates

# ==============================
# LOOP DETECTION + EXPORT
# ==============================

loop_records = []
bin_size = 2048

for idx, row in df.iterrows():
    arr = np.load(row["path"])
    arr_norm = distance_normalize(arr)

    loops = detect_loops(arr_norm, neighborhood=3, z_thresh=3.0)

    for (i, j, score) in loops:
        # remove symmetric duplicates
        if i >= j:
            continue

        g1 = row["start"] + i * bin_size
        g2 = row["start"] + j * bin_size

        start1, start2 = sorted([g1, g2])

        loop_records.append([
            row["chromosome"], start1, start1 + bin_size,
            row["chromosome"], start2, start2 + bin_size,
            score
        ])

loops_df = pd.DataFrame(loop_records, columns=[
    "Chrom1","Start1","End1","Chrom2","Start2","End2","Score"
])

# distance filtering
loops_df["distance"] = abs(loops_df["Start2"] - loops_df["Start1"])
loops_df = loops_df[loops_df["distance"] > 20000]

print("Total loops after filtering:", len(loops_df))

 ## Total loops after filtering: 12182

# export
loops_df.to_csv(bedpe_file, sep="\t", index=False, header=False)

#####


# (1) load BEDPE 

loops_df = pd.read_csv(
    bedpe_file,
    sep="\t",
    header=None,
    comment="#",
    engine="python"
)

# DEBUG: check raw structure
print("Raw shape:", loops_df.shape)
print(loops_df.head())

# Keep only rows with at least 7 columns
loops_df = loops_df.dropna(axis=0, thresh=7)

# Take first 7 columns ONLY
loops_df = loops_df.iloc[:, :7]

# Assign names
loops_df.columns = ["Chrom1","Start1","End1","Chrom2","Start2","End2","Score"]

# Convert numeric columns
for col in ["Start1","End1","Start2","End2","Score"]:
    loops_df[col] = pd.to_numeric(loops_df[col], errors="coerce")

# Drop rows where conversion failed
loops_df = loops_df.dropna()

# FINAL sanity check
print("After cleaning:", loops_df.shape)
print(loops_df.dtypes)
print(loops_df.head())


loops_df["distance"] = abs(loops_df["Start2"] - loops_df["Start1"])

print(loops_df["distance"].describe())

# keep only meaningful loops (>20 kb)
loops_df = loops_df[loops_df["distance"] > 20000]



# (2) Create PyRanges for each anchor

# Build anchor1
anchor1_df = pd.DataFrame({
    "Chromosome": loops_df["Chrom1"].astype(str),
    "Start": loops_df["Start1"],
    "End": loops_df["End1"],
    "loop_id": loops_df.index
})

# Build anchor2
anchor2_df = pd.DataFrame({
    "Chromosome": loops_df["Chrom2"].astype(str),
    "Start": loops_df["Start2"],
    "End": loops_df["End2"],
    "loop_id": loops_df.index
})

# Drop bad rows (optional but safe)
anchor1_df = anchor1_df.dropna()
anchor2_df = anchor2_df.dropna()

# Convert to PyRanges
anchor1 = pr.PyRanges(anchor1_df)
anchor2 = pr.PyRanges(anchor2_df)

# (3) Load GeneOfInterest peaks
GeneOfInterest = pr.read_bed(chip_file)

# (4) Find overlaps
a1_overlap = anchor1.join(GeneOfInterest, slack=5000)
a2_overlap = anchor2.join(GeneOfInterest, slack=5000)

# (5) Map back to loops
loops_df["anchor1_GeneOfInterest"] = loops_df.index.isin(a1_overlap.loop_id)
loops_df["anchor2_GeneOfInterest"] = loops_df.index.isin(a2_overlap.loop_id)

loops_df["any_GeneOfInterest"] = loops_df["anchor1_GeneOfInterest"] | loops_df["anchor2_GeneOfInterest"]
loops_df["both_GeneOfInterest"] = loops_df["anchor1_GeneOfInterest"] & loops_df["anchor2_GeneOfInterest"]


#####

print("Total loops:", len(loops_df))
print("GeneOfInterest loops:", loops_df["any_GeneOfInterest"].sum())
print("Fraction GeneOfInterest:", loops_df["any_GeneOfInterest"].mean())

## Total loops: 12182
## GeneOfInterest loops: 852
## Fraction GeneOfInterest: 0.06993925463799047

## meaning : Fraction of GeneOfInterest:  6.9%
##           ~7% of predicted chromatin loops involve GeneOfInterest-bound regions, 
##           indicating selective but non-global association with 3D genome structure.



## Q: Do GeneOfInterest-associated loops have higher strength?

loops_df.groupby("any_GeneOfInterest")["Score"].mean()

    # Loop strength difference
    
#    any_GeneOfInterest
#    False    7.937860
#    True     8.557711
#    Name: Score, dtype: float64
    
    # what it means: 

## Q: What fraction of loops involve GeneOfInterest?

loops_df["any_GeneOfInterest"].mean()


    # 0.06993925463799047
    # Only ~7% of loops involve GeneOfInterest
    # what it means: GeneOfInterest is not a global loop organizer

## Q: Are there GeneOfInterest–GeneOfInterest loops?”

loops_df["both_GeneOfInterest"].sum()

    # 14
    
    # meaning: How many GeneOfInterest-bound loci are interacting
    #          (potential regulatory hub)


## Is GeneOfInterest enriched at loops?
## (Are GeneOfInterest-bound regions preferentially located at loop anchors?”)

bin_size = 2048
n_iter = 1000  ## number of permutations

# ------------------------------
# Real GeneOfInterest fraction in loops
# ------------------------------
real_rate = loops_df["any_GeneOfInterest"].mean()

# ------------------------------
# Permutation test
# ------------------------------
rand_rates = []

for _ in range(n_iter):

    # random shift within genomic window
    shift = np.random.randint(-50000, 50000, len(loops_df))

    rand_start = loops_df["Start1"].values + shift

    # prevent invalid coordinates
    rand_start = np.clip(rand_start, 1, None)

    rand_df = pd.DataFrame({
        "Chromosome": loops_df["Chrom1"].astype(str),
        "Start": rand_start,
        "End": rand_start + bin_size,
        "loop_id": loops_df.index
    })

    rand_pr = pr.PyRanges(rand_df)

    # GeneOfInterest overlap
    rand_overlap = rand_pr.join(GeneOfInterest, slack=5000)

    rand_flag = loops_df.index.isin(rand_overlap.loop_id)

    rand_rates.append(rand_flag.mean())

# ------------------------------
# Summary statistics
# ------------------------------
rand_rates = np.array(rand_rates)

print("========== GeneOfInterest ENRICHMENT ==========")
print("Real GeneOfInterest fraction:", real_rate)
print("Random mean:", rand_rates.mean())
print("Random std:", rand_rates.std())

# Z-score enrichment
z = (real_rate - rand_rates.mean()) / rand_rates.std(ddof=1)

print("Z-score:", z)

# Optional empirical p-value
p_empirical = np.mean(rand_rates >= real_rate)
print("Empirical p-value:", p_empirical)
print("======================================")


# ========== GeneOfInterest ENRICHMENT ==========
# Real GeneOfInterest fraction: 0.06993925463799047
# Random mean: 0.03426013790838942
# Random std: 0.001520823456050657
# Z-score: 23.448660373586108
# Empirical p-value: 0.0
# ======================================

##  Z-score:
##  The observed GeneOfInterest enrichment is ~23 standard deviations above random expectation.

## extremely strong spatial preference
# not explainable by random positioning
# not sensitive to sampling depth 
# (tested with n_iter = 100 as well and results are consistent)

# p-value
# Even with 1000 iterations: 0 / 1000 permutations exceeded observed value
# p < 0.001

## Q: Are GeneOfInterest-bound regions preferentially located at loop anchors?
## A: Yes, strongly.

## GeneOfInterest binding sites are significantly enriched at 
# predicted loop anchors compared to locally shuffled 
# genomic positions (p < 0.001, permutation test, n = 1000), 
# suggesting preferential localization within 3D regulatory interaction hubs.

# GeneOfInterest is highly enriched at 3D regulatory interaction points
# GeneOfInterest sits in spatially active chromatin hubs
# GeneOfInterest is structurally associated with loop architecture




### Run Mann-Whitney U test

# (usage: determine if there is a significant difference between two independent samples, 
#         especially when data is not normally distributed.)

## Q: Is the observed strength difference between observed and random real?

u, p = mannwhitneyu(
    loops_df[loops_df["any_GeneOfInterest"]]["Score"],
    loops_df[~loops_df["any_GeneOfInterest"]]["Score"]
)

print("p-value:", p)

  ## 0.02352185582288974
  ## statistically significant p < 0.05
  
  ## GeneOfInterest-associated loops are significantly stronger than non-GeneOfInterest loops, 
  ## despite representing only a minority (~7%) of all loops.
  
  
  ####  figures ###
  
vals = loops_df.groupby("any_GeneOfInterest")["Score"].mean()

plt.bar(["No GeneOfInterest","GeneOfInterest"], vals)
plt.ylabel("Mean loop strength")
plt.title("GeneOfInterest-associated loops are stronger")
plt.show()
  

print("Fraction GeneOfInterest loops:", loops_df["any_GeneOfInterest"].mean())
  
