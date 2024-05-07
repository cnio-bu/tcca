import scanpy as sc
import pandas as pd
import scvi
import torch
import pymde 

scvi.settings.seed = 120394
torch.set_float32_matmul_precision("high")

## load sc data
adata = sc.read_h5ad("tcca_annotated.h5ad")

scanvi_model = scvi.model.SCANVI.load("raw_scanvi_model_hgv", adata)

## retrieve latent space
SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)

## LABEL INFERENCE and PYMDE UMAP (GPU acc.)
SCANVI_MDE_KEY = "X_scANVI_MDE"
adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCANVI_LATENT_KEY])

#  Show graph
sc.pl.embedding(
    adata, basis=SCANVI_MDE_KEY, color=["cell_type_main"], ncols=1, frameon=False
)

## FULL clustering (TME + malignants)
sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)

## write back model
#adata.write("tcca_annotated.h5ad")

## Subset malignant cells but keep latent representation for reclustering

## Odd, cannot subset if .raw is set
adata.raw = None
malignancy = adata[adata.obs["malignancy"] == True]

## redo the label inference and pymde UMAP
SCANVI_MDE_KEY = "X_scANVI_MDE"
malignancy.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(malignancy.obsm[SCANVI_LATENT_KEY])

## Perform leiden subclustering on the Latent space again
sc.pp.neighbors(malignancy, use_rep=SCANVI_LATENT_KEY)
sc.tl.leiden(malignancy, key_added="expr_clusters", resolution=1)

malignancy.write("tcca_malignant_clusters.h5ad")
