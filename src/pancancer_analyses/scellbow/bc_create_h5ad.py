import scanpy as sc
import pandas as pd
import anndata as ad
import os

os.chdir("/home/lmgonzalezb/Documents/bc-meta/SCellBow/GBM/source/")

## TARGET DATASET
bc_df = pd.read_csv("gbm_mtx_source.tsv", sep = "\t")

# Create an AnnData object
adata = ad.AnnData(X = bc_df)

adata.var.index = bc_df.index
adata.obs.index = bc_df.columns


# Load metadata
metadata = pd.read_csv("gbm_metadata_source.tsv", sep="\t")
adata.obs = metadata

# Set all drugs as highly variable
adata.var["highly_variable"] = True

# Save AnnData object as an h5ad file
adata.write("gbm_bcs_source.h5ad")


## SOURCE DATASET
adata = ad.read_h5ad("gbm_bcs_target.h5ad")
metadata = pd.read_csv("gbm_metadata_target.tsv", sep="\t")
adata.obs = metadata

# Set all drugs as highly variable
adata.var["highly_variable"] = True

# Save Anndata object as an h5ad file
adata.write("gbm_bcs_target.h5ad")
