import scanpy as sc
import pandas as pd
import anndata as ad
import os
import hotspot
import matplotlib.pyplot as plt
import pickle
import sys

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
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Create column sample_study.
adata.obs["sample_study"] = adata.obs["sample"].astype(str) + adata.obs["study"].astype(
    str
)
adata.obs["data_pmid"] = adata.obs["data_pmid"].astype(str)

# Load the anndata object witht the python embeddings
adata_int = sc.read_h5ad("tcca_annotated_clustered.h5ad")
adata.obsm["X_scANVI"] = adata_int.obsm["X_scANVI"]
adata_malignant = adata[adata.obs["malignancy"] == True]
del adata, adata_int

# Compute variances for all genes
sc.pp.highly_variable_genes(
    adata_malignant,
    flavor="seurat_v3",
    n_top_genes=10000,
    layer="counts",
    batch_key="sample_study",
    subset=True,
    #   span=0.5
)

# Create the Hotspot object and the neighborhood graph from the integrated embeddings
adata_malignant.layers["counts_csc"] = adata_malignant.layers["counts"].tocsc()
# adata_malignant.write("/storage/scratch01/users/mgonzalezb/bc-meta/integration/tcca_hotspot_in.h5ad")

# Create the hotspot object
hs = hotspot.Hotspot(
    adata_malignant,
    layer_key="counts_csc",
    model="danb",
    latent_obsm_key="X_scANVI",
    umi_counts_obs_key="nCount_RNA",
)
print("Hotspot object created")

# Remove anndata objects no longer needed
del adata_malignant

# Compute neighbourhood
hs.create_knn_graph(
    weighted_graph=False,
    n_neighbors=100,
)

print("Neighbourhood of hotspot object computed")

# Compute gene autocorrelations
hs_results = hs.compute_autocorrelations(jobs=15)
print("Gene autocorrelations computed")
print(hs_results.head(15))

# Select the genes with significant lineage autocorrelation
hs_genes = (
    hs_results.loc[hs_results.FDR < 0.05]
    .sort_values("Z", ascending=False)
    .head(500)
    .index
)

# Compute pair-wise local correlations between these genes
lcz = hs.compute_local_correlations(hs_genes, jobs=15)
print("Gene local correlations computed")


# Create modules by agglomerative clustering
modules = hs.create_modules(min_gene_threshold=15, core_only=True, fdr_threshold=0.05)

print(modules.value_counts())


os.chdir("/storage/scratch01/users/mgonzalezb/bc-meta/hotspot/")

# Plot the results of hs.create_modules
hs.plot_local_correlations(vmin=-15, vmax=15)
plt.savefig("hotspot_modules_plot.png", dpi=600, bbox_inches="tight")
plt.close()

# Store the results of the autocorrelation of genes from given module
results = hs.results.join(hs.modules)
results.to_csv("gene_modules.tsv", sep="\t")

# Hotspot can compute aggregate module scores for each cell
module_scores = hs.calculate_module_scores()
print(module_scores.head())
module_scores.to_csv("module_scores.tsv", sep="\t")