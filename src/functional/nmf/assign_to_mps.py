import scanpy as sc
import pandas as pd
import scvi
import torch
import pymde
import matplotlib.pyplot as plt
import os
import seaborn as sns

scvi.settings.seed = 120394
torch.set_float32_matmul_precision("high")

# Set working directory
os.chdir("/storage/scratch01/shared/projects/bc-meta/")

# Subset malignant cells
adata = sc.read_h5ad("single_cell/seurat/tcca/tcca_annotated_clustered.h5ad")
malignants = adata[adata.obs["malignancy"] == True]

# Load UCell scores of metaprograms(MPs) in malignant population
ucell_scores_mps = pd.read_csv(
    "functional_nmf/sample_wise/mps_ucell_scores.tsv", sep=" "
)
print(ucell_scores_mps.index.tolist() == malignants.obs.cell.tolist())
malignants.obs_names = [f"c{i}" for i in range(1, malignants.n_obs + 1)]
ucell_scores_mps.index = [f"c{i}" for i in range(1, len(ucell_scores_mps) + 1)]

# Add UCells scores of cancer state MPs (exclude lineage-specific) as a new embedding
mps = [f"MP{i}_UCell" for i in range(1, 44)]
malignants.obs = pd.concat([malignants.obs, ucell_scores_mps], axis=1)

# Plot MPs on integrated UMAP
sc.pl.embedding(
    malignants,
    basis="X_scANVI_MDE",
    color=mps[30:44],
    cmap="RdYlBu_r",
    legend_loc="on data",
    vmin=0,
    vmax=0.7,
    ncols=5,
    legend_fontsize=6,
    frameon=False,
)

plt.savefig(
    "functional_nmf/sample_wise/figures/integrated_umap_mp31-43.png",
    dpi=300,
    bbox_inches="tight",
)

# Funcional annotation of cells based on most heterogeneous metaprograms
threshold = 0.3

# Count how many metaprograms are above threshold for each cell
malignants.obs["MP_annotation"] = "Other"

# Assign "MP1" to all cells with MP1 > threshold
malignants.obs.loc[malignants.obs["MP1_UCell"] > threshold, "MP_annotation"] = "MP1"


# Assign "MP3" to all cells with MP3 > threshold
malignants.obs["MP_annotation"] = malignants.obs["MP_annotation"].astype(str)
threshold = 0.4
malignants.obs.loc[
    (malignants.obs["MP3_UCell"] > threshold)
    & (malignants.obs["MP_annotation"] == "Other"),
    "MP_annotation",
] = "MP3"

# Assign "MP8" to all cells with MP8 > threshold
malignants.obs["MP_annotation"] = malignants.obs["MP_annotation"].astype(str)
threshold = 0.4
malignants.obs.loc[
    (malignants.obs["MP8_UCell"] > threshold)
    & (malignants.obs["MP_annotation"] == "Other"),
    "MP_annotation",
] = "MP8"

# Assign "MP14" to all cells with MP8 > threshold
malignants.obs["MP_annotation"] = malignants.obs["MP_annotation"].astype(str)
threshold = 0.3
malignants.obs.loc[
    (malignants.obs["MP14_UCell"] > threshold)
    & (malignants.obs["MP_annotation"] == "Other"),
    "MP_annotation",
] = "MP14"

# Assign "MP17" to all cells with MP17 > threshold
malignants.obs["MP_annotation"] = malignants.obs["MP_annotation"].astype(str)
threshold = 0.3
malignants.obs.loc[
    (malignants.obs["MP17_UCell"] > threshold)
    & (malignants.obs["MP_annotation"] == "Other"),
    "MP_annotation",
] = "MP17"

# Assign "MP17" to all cells with MP17 > threshold
malignants.obs["MP_annotation"] = malignants.obs["MP_annotation"].astype(str)
threshold = 0.3
malignants.obs.loc[
    (malignants.obs["MP19_UCell"] > threshold)
    & (malignants.obs["MP_annotation"] == "Other"),
    "MP_annotation",
] = "MP19"

# Assign "MP27" to all cells with MP17 > threshold
malignants.obs["MP_annotation"] = malignants.obs["MP_annotation"].astype(str)
threshold = 0.3
malignants.obs.loc[
    (malignants.obs["MP27_UCell"] > threshold)
    & (malignants.obs["MP_annotation"] == "Other"),
    "MP_annotation",
] = "MP27"

# Assign "MP29" to all cells with MP17 > threshold
malignants.obs["MP_annotation"] = malignants.obs["MP_annotation"].astype(str)
threshold = 0.3
malignants.obs.loc[
    (malignants.obs["MP29_UCell"] > threshold)
    & (malignants.obs["MP_annotation"] == "Other"),
    "MP_annotation",
] = "MP29"


sc.pl.embedding(
    malignants,
    basis="X_scANVI_MDE",
    color="MP_annotation",
    legend_loc="right margin",
    size=1,
    legend_fontsize=6,
    frameon=False,
)

plt.savefig(
    "functional_nmf/sample_wise/figures/integrated_umap_maxmp.png",
    dpi=300,
    bbox_inches="tight",
)
