import scanpy as sc
import pandas as pd
import scvi
import torch
import pymde
import matplotlib.pyplot as plt
import os
import seaborn as sns
import mplscience


scvi.settings.seed = 120394
torch.set_float32_matmul_precision("high")

# Set working directory.
os.chdir("/storage/scratch01/shared/projects/bc-meta/single_cell")

# ## load sc data
adata = sc.read_h5ad("./seurat/tcca/tcca_annotated_clustered.h5ad")

adata.obs["cell_type_broad"] = adata.obs["cell_type_broad"].cat.add_categories(
    ["Unknown"]
)
not_na_levels = adata.obs["cell_type_broad"]
not_na_levels[not_na_levels.isnull()] = "Unknown"
adata.obs["cell_type_broad"] = not_na_levels

scanvi_model = scvi.model.SCANVI.load("./seurat/tcca/raw_scanvi_model", adata)
scvi_model = scvi.model.SCVI.load("./seurat/tcca/raw_scvi_model", adata)

# retrieve latent space for scvi model
SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()

# ## retrieve latent space
SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)

# ## LABEL INFERENCE and PYMDE UMAP (GPU acc.)
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
    "Unknown": "#808080",
}

sc.pl.embedding(
    adata,
    basis=SCANVI_MDE_KEY,
    color=["cell_type_broad"],
    palette=cell_type_palette,
    ncols=1,
    legend_loc="on data",
    legend_fontsize=4,
    frameon=False,
)
plt.savefig("plots/integrated_celltypes.png", dpi=600, bbox_inches="tight")

# # Color by study
# sc.pl.embedding(
#    adata, basis = SCANVI_MDE_KEY, color = ["study"], ncols = 1, frameon = False
# )
# plt.savefig("plots/integrated_study.png", dpi = 600, bbox_inches = 'tight')

# # Color by sequencing technology
# sc.pl.embedding(
#    adata, basis = SCANVI_MDE_KEY, color = ["sequencing_tech"], ncols = 1, frameon = False
# )

# plt.savefig("plots/integrated_tech.png", dpi = 600, bbox_inches = 'tight')

# # Color by genome assembly
# sc.pl.embedding(
#    adata, basis = SCANVI_MDE_KEY, color = ["genome_assembly"], ncols = 1, frameon = False
# )
# plt.savefig("plots/integrated_refgenome.png", dpi = 600, bbox_inches = 'tight')


adata = sc.read_h5ad("seurat/tcca/tcca_annotated_clustered.h5ad")

## FULL clustering (TME + malignants)
sc.pp.neighbors(adata, n_neighbors=300, use_rep=SCANVI_LATENT_KEY)

# sc.tl.leiden(adata, flavor = "igraph", n_iterations = 2, directed = False, key_added = "full_clusters")

## Malignant clusters
resolutions = [0.3, 0.5, 1.0, 1.5]
for res in resolutions:
    key = f"malignant_res{res}"
    sc.tl.leiden(
        adata,
        resolution=res,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        restrict_to=("cell_type_broad", ["Malignant"]),
        key_added=key,
    )
## write back model
adata.write("tcca_annotated_clustered.h5ad")

# Plot only malignant cell clusters
adata.raw = None
malignancy = adata[adata.obs["malignancy"] == True]

# Remove clusters with only one cell
cluster_size = malignancy.obs.malignant_clusters.value_counts()
valid_clusters = cluster_size[cluster_size >= 100].index
malignancy = malignancy[malignancy.obs["malignant_clusters"].isin(valid_clusters)]

unique_clusters = pd.Categorical(malignancy.obs.malignant_clusters).categories
cluster_mapping = {old: new for new, old in enumerate(unique_clusters, start=1)}
malignancy.obs["malignant_clusters"] = malignancy.obs["malignant_clusters"].map(
    cluster_mapping
)

# Plot malignant clusters
sc.pl.embedding(
    malignancy,
    basis=SCANVI_MDE_KEY,
    color=["malignant_clusters"],
    legend_loc="on data",
    ncols=1,
    legend_fontsize=6,
    frameon=False,
)
plt.savefig("plots/integrated_malignant_clusters.png", dpi=600, bbox_inches="tight")


# Check annotation using cell type markers:
marker_genes_dict = [
    "CD19",
    "CD79A",
    "MS4A1",
    "IGKC",
    "IGHG1",
    "CD4",
    "ILR7",
    "TRAC",
    "CD3D",
    "CD8A",
    "CD8B",
    "CCL5",
    "LAG3",
    "FOXP3",
    "CTLA4",
    "CD127",
    "GATA3",
    "SOX4",
    "TRDV1",
    "NKX1-1",
    "CD74",
    "HLA-DRA",
    "CD11c",
    "CLEC4C",
    "IRF7",
    "IL3RA",
    "AZU1",
    "CEACAM8",
    "CTSS",
    "NOS2",
    "LYZ",
    "CD14",
    "CD68",
    "NKG7",
    "KLRD1",
    "GZMB",
    "CPA3",
    "TPSAB1",
    "TPSB2",
    "HBA1",
    "AHSP",
    "PPBP",
    "PF4",
    "ITGA2B",
    "EPCAM",
    "KRT8",
    "MUC1",
    "PECAM1",
    "ACKR1",
    "ABCG2",
    "COL1A1",
    "VIM",
    "PDGFRA",
    "CD44",
    "CD34",
    "ENG",
    "GFAP",
    "OLIG2",
    "PLP1",
    "MAP2",
    "RBFOX3",
    "SYN1",
]

marker_genes_dict = {
    "B-cell": ["CD19", "CD79A", "MS4A1"],
    "Plasma cell": ["IGKC", "IGHG1"],
    "CD4+ T-cell": ["CD4", "ILR7", "TRAC", "CD3D"],
    "CD8+ T-cell": ["CD8A", "CD8B", "CCL5"],
    "Regulatory T-cell": ["LAG3", "FOXP3", "CTLA4", "NKX1-1"],
    "Innate lymphoid cells": ["CD127", "GATA3", "SOX4"],
    "Unconventional T-cells": ["KLRK1", "TRDV1", "NKX1-1"],
    "NK cell": ["NKG7", "KLRD1", "GZMB"],
    "Dendritic cell": ["CD74", "HLA-DRA", "CD11c"],
    "Plasmacytoid dendritic cell": ["CLEC4C", "IRF7", "IL3RA"],
    "Granulocyte": ["AZU1", "CEACAM8", "CTSS"],
    "Monocyte/Macrophage": ["CTSS", "NOS2", "LYZ", "CD14", "CD68"],
    "Mast": ["CPA3", "TPSAB1", "TPSB2"],
    "Erythrocyte": ["HBA1", "AHSP"],
    "Platelet": ["PPBP", "PF4", "ITGA2B"],
    "Epithelial": ["EPCAM", "KRT8", "MUC1"],
    "Endothelial": ["PECAM1", "ACKR1", "ABCG2"],
    "Stromal cell": ["COL1A1", "VIM", "PDGFRA"],
    "Stem": ["CD44", "CD34", "ENG"],
    "Glial cell": ["GFAP", "OLIG2", "PLP1"],
    "Neuron": ["MAP2", "RBFOX3", "SYN1"],
}

filtered_marker_genes_dict = {}

# Iterate through the original dictionary
for cell_type, genes in marker_genes_dict.items():
    # Keep only genes that are in the valid_genes list
    filtered_genes = [gene for gene in genes if gene in adata.var_names]

    # Only add to the new dictionary if there are valid genes left
    if filtered_genes:
        filtered_marker_genes_dict[cell_type] = filtered_genes

sc.pl.stacked_violin(
    adata,
    filtered_marker_genes_dict,
    groupby="cell_type_broad",
    swap_axes=False,
    dendrogram=True,
)

sc.pl.matrixplot(
    adata,
    filtered_marker_genes_dict,
    groupby="cell_type_broad",
    categories_order=filtered_marker_genes_dict.keys(),
    dendrogram=False,
    cmap="Blues",
    standard_scale="var",
    colorbar_title="column scaled\nexpression",
)

plt.savefig("plots/marker_genes_celltype.png", dpi=600, bbox_inches="tight")


# # Plot Hotspot gene modules
# adata_malignant = adata[adata.obs["malignancy"] == True]
# module_scores = pd.read_csv("hotspot/module_scores.tsv", sep = "\t", index_col = "cell_id")
# module_scores.index = module_scores.index.astype(str)
# adata_malignant.obs = pd.merge(adata_malignant.obs,
#                                module_scores,
#                                how = "inner",
#                                left_index = True,
#                                right_index = True)

# plt.figure(figsize=(10, 10))
# with mplscience.style_context():
#     sc.pl.embedding(adata_malignant, basis = SCANVI_MDE_KEY, color = module_scores.columns,
#     frameon=False, vmin=-1, vmax=1
# )
# plt.show()
# plt.savefig("hotspot/gene_modules_umap.png", dpi = 600, bbox_inches = 'tight')

# # Subset malignant cells but keep latent representation for reclustering


# ## Odd, cannot subset if .raw is set
# # adata.raw = None
# # malignancy = adata[adata.obs["malignancy"] == True]

# ## redo the label inference and pymde UMAP
# #SCANVI_MDE_KEY = "X_scANVI_MDE"
# #malignancy.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(malignancy.obsm[SCANVI_LATENT_KEY])

# ## Perform leiden subclustering on the Latent space again
# #sc.pp.neighbors(malignancy, use_rep=SCANVI_LATENT_KEY)
# #sc.tl.leiden(malignancy, key_added="expr_clusters", resolution=1)

# #malignancy.write("tcca_malignant_clusters.h5ad")
