import scanpy as sc
import pandas as pd
import scvi
import torch

scvi.settings.seed = 120394
print("Last run with scvi-tools version:", scvi.__version__)

torch.set_float32_matmul_precision("high")

adata = sc.read_h5ad("tcca_raw_mat.h5ad")
print("File loaded!")

print("Load metadata")
metadata = pd.read_csv("tcca_annotation_raw.tsv", sep="\t")
metadata.cell_id = metadata.cell_id.astype(str)
metadata.set_index("cell_id", inplace=True)

print("Metadata set!")
adata.obs = metadata


scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=["sample", "study"])
print("Generated model load")

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()
print("Done! saving raw model")
model.save("./raw_scvi_model", overwrite=True)

print("Generate SCANVI model")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=adata,
    labels_key="cell_type_main",
    unlabeled_category="unknown",
)
scanvi_model.train(max_epochs=20, n_samples_per_label=100)
print("DONE! saving...")
scanvi_model.save("./raw_scanvi_model", overwrite=True)