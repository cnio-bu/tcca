import scanpy as sc
import pandas as pd
import scvi
import torch
import anndata as ad
import os

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
adata.obs["sample_study"] = adata.obs["sample"].astype(str) + adata.obs["study"].astype(str)
adata.obs["data_pmid"] = adata.obs["data_pmid"].astype(str)

# Keep full dimension anndata object safe.
adata.raw = adata

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key="sample_study",
    subset=True,
    #   span=0.5
)
adata.write(
    "/storage/scratch01/users/mgonzalezb/bc-meta/integration/tcca_annotated.h5ad"
)

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample_study")

print("Generated model load")

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()
print("Done! saving raw model")
model.save("./raw_scvi_model", overwrite=True)

print("Generate SCANVI model")

## TODO: make sure no "NaN" in main cell type. Set to
## "unknown" if needed
adata.obs["cell_type_main"] = adata.obs["cell_type_main"].cat.add_categories(
    ["unknown"]
)
not_na_levels = adata.obs["cell_type_main"]
not_na_levels[not_na_levels.isnull()] = "unknown"
adata.obs["cell_type_main"] = not_na_levels

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=adata,
    labels_key="cell_type_main",
    unlabeled_category="unknown",
)
scanvi_model.train(max_epochs=20, n_samples_per_label=100)
print("DONE! saving...")
scanvi_model.save("./raw_scanvi_model", overwrite=True)
