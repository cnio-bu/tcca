import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
from scib_metrics.nearest_neighbors import NeighborsResults
import faiss
import time

# Set working directory.
os.chdir("/storage/scratch01/shared/projects/bc-meta/single_cell")

# Load sc data
adata = sc.read_h5ad("./seurat/tcca/tcca_annotated_clustered.h5ad")

sc.pp.subsample(adata, n_obs=500_000)

# Compute pca
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
plt.savefig("./seurat/tcca/plots/pca_variance_ratio.png", 
            dpi = 600, 
            bbox_inches = 'tight')
adata.obsm["Unintegrated"] = adata.obsm["X_pca"]

# Custom nearest neighbour function from the scib-metrics tutorial to speep up
# the computation of the metrics
def faiss_hnsw_nn(X: np.ndarray, k: int):
    """Gpu HNSW nearest neighbor search using faiss.

    See https://github.com/nmslib/hnswlib/blob/master/ALGO_PARAMS.md
    for index param details.
    """
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    M = 32
    index = faiss.IndexHNSWFlat(X.shape[1], M, faiss.METRIC_L2)
    gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    gpu_index.add(X)
    distances, indices = gpu_index.search(X, k)
    del index
    del gpu_index
    # distances are squared
    return NeighborsResults(indices=indices, distances=np.sqrt(distances))


def faiss_brute_force_nn(X: np.ndarray, k: int):
    """Gpu brute force nearest neighbor search using faiss."""
    X = np.ascontiguousarray(X, dtype=np.float32)
    res = faiss.StandardGpuResources()
    index = faiss.IndexFlatL2(X.shape[1])
    gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    gpu_index.add(X)
    distances, indices = gpu_index.search(X, k)
    del index
    del gpu_index
    # distances are squared
    return NeighborsResults(indices=indices, distances=np.sqrt(distances))


# Benchmarking integration methods
biocons = BioConservation(isolated_labels=False)
start = time.time()
bm = Benchmarker(
    adata,
    batch_key="sample_study",
    label_key="cell_type_main",
    embedding_obsm_keys=["Unintegrated", "X_scANVI", "X_scVI"],
    pre_integrated_embedding_obsm_key="X_pca",
    bio_conservation_metrics=biocons,
    batch_correction_metrics=BatchCorrection(),
    n_jobs=-1,
)
bm.prepare(neighbor_computer=faiss_brute_force_nn)
bm.benchmark()
end = time.time()
print(f"Time: {int((end - start) / 60)} min {int((end - start) % 60)} sec")

bm.plot_results_table()
plt.savefig("./seurat/tcca/plots/integration_benchmark_scaled.png", 
            dpi = 600, 
            bbox_inches = 'tight')

bm.plot_results_table(min_max_scale=False)
plt.savefig("./seurat/tcca/plots/integration_benchmark.png", 
            dpi = 600, 
            bbox_inches = 'tight')

df = bm.get_results(min_max_scale=False)
df.to_csv("./seurat/tcca/integration_benchmark.tsv", sep='\t')