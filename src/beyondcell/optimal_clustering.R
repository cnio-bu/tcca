library(tidyverse)
library(BPCells)
library(Seurat)
library(clustree)
library(bluster)
library(clv)
library(future.apply)
library(Rfast)
library(parallel)

## Set options
options(future.globals.maxSize = 20 * 1024 ^ 3)
options(Seurat.object.assay.version = 'v5')

setwd("/storage/scratch01/shared/projects/bc-meta/beyondcell")

seu <- readRDS("beyondcell_pancancer.rds")
seu <- FindVariableFeatures(object = seu,
                            selection.method = "disp",
                            nfeatures = 300)

seu <- ScaleData(seu,
                 do.scale = TRUE,
                 do.center = TRUE,
                 scale.max = 20)

seu <- RunPCA(seu, reduction.name = "pca", npcs = 100)

seu <- seu[, rownames(seu@reductions$pca@cell.embeddings)]


# Compute distance matrix for silhoutte score.
dist.matrix <- Rfast::Dist(seu@reductions$pca@cell.embeddings)

# Select the optimal number of neighbours (k) and resolution (res) using different methods.
evaluate_clustering <- function(data,
                                obj,
                                membership,
                                dist.matrix,
                                k.param,
                                res) {
  intraclust = c("centroid")
  interclust = c("centroid")
  
  if (length(unique(membership)) > 1) {
    clust_pred <- as.integer(membership)
    
    # compute intercluster distances and intracluster diameters
    cls.scatt <- cls.scatt.data(data, clust_pred, dist = "euclidean")
    
    # compute Davies-Bouldin index
    dbi_out <- clv.Davies.Bouldin(cls.scatt, intraclust, interclust)
    
    # compute silhoutte scores
    sil <- silhouette(x = clust_pred, dist = dist.matrix)
    sil <- mean(sil[, 3])
    
    # cluster purity
    purity <- neighborPurity(data, clust_pred)
    purity <- mean(purity[, 1])
    
    # cluster stability
    stab <- mclapply(1:20, function(i) {
      bootstrap <- sample(colnames(obj),
                          size = 0.8 * ncol(obj),
                          replace = FALSE)
      bootstrap_seu <- obj[, bootstrap]
      bootstrap_seu <- FindNeighbors(
        bootstrap_seu,
        dims = 1:30,
        reduction = "pca",
        k.param = 20
      )
      bootstrap_seu <- FindClusters(bootstrap_seu,
                                    cluster.name = "clusters",
                                    resolution = 0.1)
      
      # Calculate cluster stability using pairwiseRand
      clust_bootstrap <- bootstrap_seu@meta.data[bootstrap, "clusters"]
      stab <- pairwiseRand(as.factor(clust_pred[bootstrap]), clust_bootstrap, mode = "index")
      return(stab)
    }, mc.cores = 20)
    stab <-  mean(na.omit(stab))
    print(stab)
    return(list(
      dbi = dbi_out,
      sil = sil,
      purity = purity,
      stab = stab
    ))
  } else{
    return(NULL)
  }
}


metrics_tc <- lapply(c(20, 30, 50, 75, 100, 150, 200, 300), function(X) {
  print(paste0("k neighbours = ", X))
  obj <- FindNeighbors(seu,
                       dims = 1:50,
                       reduction = "pca",
                       k.param = X)
  metrics_res <- lapply(seq(.1, 1, 0.1), function(Y) {
    obj <- FindClusters(obj, cluster.name = "clusters", resolution = Y)
    metrics <- evaluate_clustering(
      data = seu@reductions$pca@cell.embeddings,
      obj = seu,
      membership = obj$seurat_clusters,
      dist.matrix = dist.matrix,
      k.param = X,
      res = Y
    )
    print(Y)
    return(metrics)
  })
  return(metrics_res)
})


# Save the results
saveRDS(metrics_tc, "clustering_metrics.rds")

# Set the names based on k neighbour and resolution parameter
names(metrics_tc) <- paste0("k.", c(20, 30, 50, 75, 100, 150, 200, 300))
metrics_tc <- lapply(metrics_tc, function(X) {
  names(X) <- paste0("res.", seq(.1, 1, 0.1))
  return(X)
})

# Extract parameter values with lower Davies Bouldin Index.
dbi <- lapply(metrics_tc, function(X) {
  dbi <- lapply(X, function(Y) {
    return(Y$dbi)
  })
  return(dbi)
})

dbi <- unlist(dbi)
k_dbi <- unlist(strsplit(names(which.min(dbi)), split = ".", fixed = TRUE))[2]
res_dbi <- paste(unlist(strsplit(
  names(which.min(dbi)), split = ".", fixed = TRUE
))[4:5], collapse = ".")


# Extract parameter values with higher silhoutte scores.
sil <- lapply(metrics_tc, function(X) {
  sil <- lapply(X, function(Y) {
    return(Y$sil)
  })
  return(sil)
})

sil <- unlist(sil)
k_sil <- unlist(strsplit(names(which.max(sil)), split = ".", fixed = TRUE))[2]
res_sil <- paste(unlist(strsplit(
  names(which.max(sil)), split = ".", fixed = TRUE
))[4:5], collapse = ".")


# Extract parameter values with higher purity scores.
purity <- lapply(metrics_tc, function(X) {
  purity <- lapply(X, function(Y) {
    return(Y$purity)
  })
  return(purity)
})

purity <- unlist(purity)
k_purity <- unlist(strsplit(names(which.max(purity)), split = ".", fixed = TRUE))[2]
res_purity <- paste(unlist(strsplit(
  names(which.max(purity)), split = ".", fixed = TRUE
))[4:5], collapse = ".")

metrics_tc <- readRDS("clustering_metrics.rds")
names(metrics_tc) <- paste0("k.", c(20, 30, 50, 75, 100, 150, 200, 300))
metrics_tc <- lapply(metrics_tc, function(X) {
  names(X) <- paste0("res.", seq(.1, 1, 0.1))
  return(X)
})

metrics_tc <- do.call(c, metrics_tc)
metrics_tc <- lapply(metrics_tc, unlist)
metrics_tc <- do.call(rbind, metrics_tc)
params <- as.data.frame(metrics_tc) %>%
  arrange(dbi, desc(sil), desc(purity)) %>%
  slice_head(n = 10) %>%
  row.names()

stability <- c()
for (x in params) {
  k.param <- unlist(strsplit(x, split = ".", fixed = TRUE))[2]
  print(k.param)
  res <- paste0(unlist(strsplit(x, split = ".", fixed = TRUE))[4:5], collapse = ".")
  print(res)
  
  stab <- mclapply(1:20, function(i) {
    bootstrap <- sample(colnames(seu),
                        size = 0.8 * ncol(seu),
                        replace = FALSE)
    bootstrap_seu <- seu[, bootstrap]
    bootstrap_seu <- FindNeighbors(
      bootstrap_seu,
      dims = 1:30,
      reduction = "pca",
      k.param = k.param
    )
    bootstrap_seu <- FindClusters(bootstrap_seu,
                                  cluster.name = "clusters",
                                  resolution = res)
    
    # Calculate cluster stability using pairwiseRand
    clust_bootstrap <- bootstrap_seu@meta.data[bootstrap, "clusters"]
    stab <- pairwiseRand(as.factor(clust_pred[bootstrap]), clust_bootstrap, mode = "index")
    return(stab)
  }, mc.cores = 20)
  stab <-  mean(unlist(stab))
  stability <- c(stability, stab)
}
