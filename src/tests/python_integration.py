import scanpy as sc
import scvi
import torch

scvi.settings.seed = 120394
print("Last run with scvi-tools version:", scvi.__version__)

torch.set_float32_matmul_precision("high")

adata = sc.read_h5ad("tcca_raw_mat.h5ad")
print("File loaded!")

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=["sample", "study"])
print("Generated model load")

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()
print("Done! saving raw model")
model.save("./raw_scvi_model", overwrite=True)