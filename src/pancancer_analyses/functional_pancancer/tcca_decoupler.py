import scanpy as sc
import decoupler as dc
import pandas as pd
import os 

# Set working directory.
os.chdir("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/")

# Load raw count matrix.
adata = sc.read_h5ad("tcca_integrated.h5ad")
print("File loaded!")

# Get the Progeny curated collection of pathways and their target genes, 
# with weights for each interaction
progeny = dc.get_progeny(organism = "human", top = 500)
progeny

print("Progeny signatures loaded")
adata.raw = adata

# Infer pathway enrichment scores with the multivariate linear model (mlm) method
dc.run_mlm(
    mat = adata,
    net = progeny,
    source = "source",
    target = "target",
    weight = "weight",
    verbose = True
)

# Visualize the obtained scores
acts = dc.get_acts(adata, obsm_key = "mlm_estimate")
acts

acts.write_h5ad("tcca_decoupler.h5ad")
