import scanpy as sc
import pandas as pd
import scvi
import torch
import anndata as ad
import scanorama
import os
import matplotlib.pyplot as plt
import numpy as np

scvi.settings.seed = 120394
print("Last run with scvi-tools version:", scvi.__version__)

torch.set_float32_matmul_precision("high")

# Set working directory.
os.chdir("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/")

# Load raw count matrix.
adata = sc.read_h5ad("tcca_raw_mat.h5ad")
print("File loaded!")

# Add the metadata from tsv file to anndata object.
print("Load metadata")
metadata = pd.read_csv("tcca_annotation_raw.tsv", sep = "\t")
metadata.cell_id = metadata.cell_id.astype(str)

metadata.set_index("cell_id", inplace = True)

print("Metadata set!")
adata.obs = metadata

# Add counts to a layer and normalize.
print("Set layer data")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)

# Create column sample_study and modify types
adata.obs["sample_study"] = adata.obs["sample"].astype(str) + adata.obs["study"].astype(str)
adata.obs["sample_study"] = adata.obs["sample_study"].astype("category")
adata.obs["data_pmid"] = adata.obs["data_pmid"].astype(str)


# Plot marker genes in a heatmap
# marker_genes_dict = {
#     "B-cell": ["CD19", "CD79A", "MS4A1"],
#     "Plasma cell": ["PRDM1", "JCHAIN", "PLA2G2D"],
#     "CD4+ T-cell": ["CD40LG", "CD3D", "IL7R"],
#     "CD8+ T-cell": ["CD8A", "CD8B", "CCL5"],
#     "Regulatory T-cell": ["FOXP3", "CTLA4"],
#     "Innate lymphoid cells": ["IL4I1", "NFKBIA", "CD69"],
#     "Unconventional T-cells": ["KLRK1", "KLRG1"],
#     "NK cell": ["NKG7", "KLRD1"],
#     "Dendritic cell": ["CD1C", "CLEC10A"],
#     "Plasmacytoid dendritic cell": ["CXCR3", "IRF7"],
#     "Granulocyte": ["NTRK1", "CTSS", "MPO", "SIGLEC8"],
#     "Monocyte/Macrophage": ["FCN1", "CSF1R", "CD68", "MSR1"],
#     "Mast": ["CPA3", "TPSAB1", "TPSB2"],
#     "Erythrocyte": ["HBA1", "HBD", "AHSP"],
#     "Platelet": ["GP9", "PF4"],
#     "Epithelial": ["EPCAM", "KRT8", "MUC1"],
#     "Endothelial": ["PECAM1", "ACKR1", "ABCG2"],
#     "Stromal cell": ["COL1A1", "DCN", "PDPN", "FAP"],
#     "Stem": ["NANOG", "ITLN2"],
#     "Glial cell": ["GFAP", "OLIG2", "PLP1"],
#     "Neuron": ["MAP2", "RBFOX3", "NXPH2", "SYNPR"],
#     "Malignant":[],
#     "Unknown":[]
# }

# sc.pl.matrixplot(
#     adata,
#     dict(list(marker_genes_dict.items())[:-2]),
#     groupby = "cell_type_broad",
#     categories_order = pd.Categorical(list(marker_genes_dict.keys())),
#     #dendrogram=True,
#     standard_scale="var",
#     #cmap = "Blues",
#     figsize=(14, 6),
#     colorbar_title="Scaled gene\nexpression",
# )


# plt.savefig("plots/marker_genes_celltype_viridis.png", dpi = 600, bbox_inches = 'tight')

# Keep full dimension anndata object safe.
adata.raw = adata

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer = "counts",
    batch_key = "sample_study",
    subset = True
 #   span=0.5
    )


scvi.model.SCVI.setup_anndata(adata, layer = "counts", batch_key = "sample_study")

print("Generated model load")

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train(max_epochs = 50)
print("Done! saving raw model")
model.save("./new_scvi_model", overwrite=True)


#  Make sure no "NaN" in main cell type. Set to "unknown" if needed
adata.obs.loc[adata.obs["cell_type_broad"].isnull(), "cell_type_broad"] = "Unknown"
adata.obs["cell_type_broad"] = adata.obs["cell_type_broad"].astype("category")

print("Generate SCANVI model")

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata = adata,
    labels_key = "cell_type_broad",
    unlabeled_category = "Unknown"
)
scanvi_model.train(max_epochs = 50, n_samples_per_label = 100)
print("DONE! saving...")
scanvi_model.save("./new_scanvi_model", overwrite = True)

# Run scanorama integration
import scanorama

# List of adata per batch
batch_cats = adata.obs.sample_study.cat.categories
adata_list = [adata[adata.obs.sample_study == b].copy() for b in batch_cats]

# Compute scanorama integration
scanorama.integrate_scanpy(adata_list)
print("Running scanorama integration")

adata.obsm["Scanorama"] = np.zeros((adata.shape[0], adata_list[0].obsm["X_scanorama"].shape[1]))
for i, b in enumerate(batch_cats):
    adata.obsm["Scanorama"][adata.obs.sample_study == b] = adata_list[i].obsm["X_scanorama"]

adata.write("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/tcca_annotated.h5ad")