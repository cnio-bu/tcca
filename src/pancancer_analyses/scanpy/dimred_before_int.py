import scanpy as sc
import pandas as pd
import anndata as ad
import os
import matplotlib.pyplot as plt


# Set working directory.
os.chdir("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/")

# Load raw count matrix.
adata = sc.read_h5ad("tcca_raw_mat.h5ad")
print("File loaded!")

# Add the metadata from tsv file to anndata object.
print("Load metadata")
metadata = pd.read_csv("tcca_annotation_raw.tsv", sep="\t")
metadata.cell_id = metadata.cell_id.astype(str)
metadata.set_index("cell_id", inplace=True)

print("Metadata set!")
adata.obs = metadata

# Add counts to a layer and normalize.
print("Set layer data")
adata.layers["counts"] = adata.X.copy()

# Make sure data is normalised, and log1p transformed before identifying highly variable genes
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Create column sample_study.
adata.obs["sample_study"] = adata.obs["sample"].astype(str) + adata.obs["study"].astype(str)
adata.obs["data_pmid"] = adata.obs["data_pmid"].astype(str)

# Select highly variable genes by normalized variance
# (seurat flavor is the default, which expects log1p normalized data)
sc.pp.highly_variable_genes(
    adata, n_top_genes=int(0.1 * adata.n_vars), flavor="seurat", subset=True
)

# Run PCA, we calculate 100 PC components here to "over-compute" for demonstration pourposes.
sc.tl.pca(adata, use_highly_variable=True, n_comps=100)

# Visualise the percent variance explained by each PC
# This is called a PCA scree plot, we can see that it reaches a plateau at around 80 PCs
sc.pl.pca_variance_ratio(adata, n_pcs=100, log=True)

# Save the figure
plt.savefig("plots/pca_variance_ratio.png", dpi=300)

# First, we need to compute the neighbour graph of cells based on their PCs.
sc.pp.neighbors(adata, n_pcs=80)

# Now, we visualize the neighbour graph by calculating UMAP.
sc.tl.umap(adata)


# Set color palette for cell types.
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
    "unknown": "#808080",
}

sc.pl.umap(adata, color="cell_type_broad", palette=cell_type_palette, frameon=False)
plt.savefig("plots/unintegrated_celltypes.png", dpi=600, bbox_inches="tight")

# Color umap by study
sc.pl.umap(adata, color="study", frameon=False)
plt.savefig("plots/unintegrated_study.png", dpi=600, bbox_inches="tight")

# Color by sequencing technology
sc.pl.umap(adata, color="sequencing_tech", frameon=False)
plt.savefig("plots/unintegrated_tech.png", dpi=600, bbox_inches="tight")


# Color umap by genome assembly
sc.pl.umap(adata, color="genome_assembly", frameon=False)
plt.savefig("plots/unintegrated_refgenome.png", dpi=600, bbox_inches="tight")


# save processed adata
adata.write_h5ad("dim_reduced_unintegrated.h5ad")
