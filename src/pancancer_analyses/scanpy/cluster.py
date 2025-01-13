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

# Set working directory.
os.chdir("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/")

## load sc data
adata = sc.read_h5ad("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/tcca_annotated.h5ad")

adata.obs["cell_type_main"] = adata.obs["cell_type_main"].cat.add_categories(["unknown"])
not_na_levels = adata.obs["cell_type_main"]
not_na_levels[not_na_levels.isnull()] = "unknown"
adata.obs["cell_type_main"] = not_na_levels

scanvi_model = scvi.model.SCANVI.load("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/raw_scanvi_model", adata)

## retrieve latent space
SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)

## LABEL INFERENCE and PYMDE UMAP (GPU acc.)
SCANVI_MDE_KEY = "X_scANVI_MDE"
adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCANVI_LATENT_KEY])

# Plot integrated umap color by cell type
cell_type_palette = {
    "Stromal cell": "#d7bafd", 
    "Glial cell": "#873e23", 
    "Epithelial": "#006cd1", 
    "B-cell": "#935ee7", 
    "Dendritic cell": "#009f83", 
    "Plasmacytoid dendritic cell": "#b3d266", 
    "Endothelial": "#d20086", 
    "Erythrocyte": "#7dda84", 
    "Granulocyte": "#ff44d4", 
    "Mast": "#0092f4", 
    "Platelet": "#df8f00", 
    "Monocyte/Macrophage": "#b66c00", 
    "NK cell": "#b88c9c", 
    "Neuron": "#cb4dd0", 
    "Stem": "#5fcaff", 
    "Innate lymphoid cells": "#7591ff", 
    "Plasma cell": "#f97d1b", 
    "CD4+ T-cell": "#394f94", 
    "CD8+ T-cell": "#d1da32", 
    "Regulatory T-cell": "#8e2c68", 
    "Unconventional T-cells": "#3d6700", 
    "Malignant": "#d82500",
     "unknown": "#808080"
}

sc.pl.embedding(
   adata, basis = SCANVI_MDE_KEY, color=["cell_type_broad"], palette = cell_type_palette, ncols = 1, frameon = False
)
plt.savefig("plots/integrated_celltypes.png", dpi = 600, bbox_inches = 'tight')

# Color by study
sc.pl.embedding(
   adata, basis = SCANVI_MDE_KEY, color = ["study"], ncols = 1, frameon = False
)
plt.savefig("plots/integrated_study.png", dpi = 600, bbox_inches = 'tight')

# Color by sequencing technology
sc.pl.embedding(
   adata, basis = SCANVI_MDE_KEY, color = ["sequencing_tech"], ncols = 1, frameon = False
)

plt.savefig("plots/integrated_tech.png", dpi = 600, bbox_inches = 'tight')

# Color by genome assembly
sc.pl.embedding(
   adata, basis = SCANVI_MDE_KEY, color = ["genome_assembly"], ncols = 1, frameon = False
)
plt.savefig("plots/integrated_refgenome.png", dpi = 600, bbox_inches = 'tight')

## FULL clustering (TME + malignants)
sc.pp.neighbors(adata, use_rep = SCANVI_LATENT_KEY)

sc.tl.leiden(adata, flavor = "igraph", n_iterations = 2, directed = False, key_added = "full_clusters")
## Malignant clusters
sc.tl.leiden(adata, flavor = "igraph", n_iterations = 2, directed = False, restrict_to = ("cell_type_main", ["malignant"]), key_added = "malignant_clusters")
## write back model
adata.write("tcca_annotated_clustered.h5ad")

## Subset malignant cells but keep latent representation for reclustering

## Odd, cannot subset if .raw is set
#adata.raw = None
#malignancy = adata[adata.obs["malignancy"] == True]

## redo the label inference and pymde UMAP
#SCANVI_MDE_KEY = "X_scANVI_MDE"
#malignancy.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(malignancy.obsm[SCANVI_LATENT_KEY])

## Perform leiden subclustering on the Latent space again
#sc.pp.neighbors(malignancy, use_rep=SCANVI_LATENT_KEY)
#sc.tl.leiden(malignancy, key_added="expr_clusters", resolution=1)

#malignancy.write("tcca_malignant_clusters.h5ad")
