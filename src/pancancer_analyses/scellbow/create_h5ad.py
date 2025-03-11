import scanpy as sc
import pandas as pd
import anndata as ad
import os

os.chdir("/storage/scratch01/shared/projects/bc-meta/single_cell/scellbow/")

############################### BEYONDCELL DATA ###############################
## Source dataset
bc_df = pd.read_csv("beyondcell/gbm_mtx_source.tsv", sep = "\t")

# Create an AnnData object
adata = ad.AnnData(X = bc_df.transpose())

adata.var.index = bc_df.index
adata.obs.index = bc_df.columns


# Load metadata
metadata = pd.read_csv("beyondcell/gbm_metadata_source.tsv", sep="\t")
adata.obs = metadata

# Set all drugs as highly variable
adata.var["highly_variable"] = True

# Save AnnData object as an h5ad file
adata.write("beyondcell/gbm_bcs_source.h5ad")


## Target dataset
adata = ad.read_h5ad("beyondcell/gbm_bcs_target.h5ad")
metadata = pd.read_csv("beyondcell/gbm_metadata_target.tsv", sep="\t")
adata.obs = metadata

# Set all drugs as highly variable
adata.var["highly_variable"] = True

# Save Anndata object as an h5ad file
adata.write("beyondcell/gbm_bcs_target.h5ad")


############################### EXPRESSION DATA ###############################
## Source dataset
expr_df = pd.read_csv("gbm_mtx_source_RNA.tsv", sep="\t")

# Create an AnnData object
adata = ad.AnnData(X=expr_df.transpose())

adata.var.index = expr_df.index
adata.obs.index = expr_df.columns

# Load metadata
metadata = pd.read_csv("gbm_metadata_source_RNA.tsv", sep="\t")
adata.obs = metadata

# Save AnnData object as an h5ad file
adata.write("expression/gbm_expr_source.h5ad")


## Target dataset
adata = ad.read_h5ad("expression/gbm_expr_target.h5ad")
metadata = pd.read_csv("expression/gbm_metadata_target.tsv", sep="\t")
adata.obs = metadata

# Save Anndata object as an h5ad file
adata.write("expression/gbm_expr_target.h5ad")


