library(dplyr)
library(tidyverse)
library(circlize)
library(factoextra)
set.seed(123)
setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/sctherapy/")
source("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")

similarity_matrix <- readRDS("results/jaccard_matrix.rds")
cluster_assignment <- readRDS("results/speclustering_reordered.rds")
dist_matrix <- as.dist(1 - similarity_matrix)

# 1. ARI spectral vs hierarchical
hc <- hclust(dist_matrix, method = "ward.D2")
hc_clusters <- cutree(hc, k = 10)
ari <- adjustedRandIndex(as.integer(cluster_assignment), hc_c∫lusters)
cat("ARI spectral vs hierarchical:", ari, "\n")

# 2. Silhouette score para k=2:14
k_range <- 2:14
silhouette_scores <- sapply(k_range, function(k) {
  cat("Running specc k =", k, "\n")
  result <- specc(similarity_matrix, centers = k)
  sil <- silhouette(as.integer(result), dist_matrix)
  mean(sil[, 3])
})

# k=10 con tu clustering ya calculado (evita recalcular)
sil_k10 <- silhouette(as.integer(cluster_assignment), dist_matrix)
cat("Silhouette k=10 (tu clustering):", mean(sil_k10[, 3]), "\n")

# Plot
sil_df <- data.frame(k = k_range, silhouette = silhouette_scores)
ggplot(sil_df, aes(x = k, y = silhouette)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = 10, col = "red", linetype = "dashed") +
  labs(x = "Number of clusters (k)", y = "Mean silhouette score",
       title = "Silhouette score across k") +
  theme_classic()
ggsave("figures/silhouette_scores.pdf", width = 6, height = 4)

# 3. Within-cluster Jaccard coherence
compute_coherence <- function(clusters, sim_matrix) {
  mean(sapply(unique(clusters), function(cl) {
    idx <- which(clusters == cl)
    if (length(idx) < 2) return(NA)
    vals <- sim_matrix[idx, idx]
    mean(vals[lower.tri(vals)])
  }), na.rm = TRUE)
}

coherence_scores <- sapply(k_range, function(k) {
  result <- specc(similarity_matrix, centers = k)
  compute_coherence(as.integer(result), similarity_matrix)
})

coherence_k10 <- compute_coherence(as.integer(cluster_assignment), similarity_matrix)
cat("Within-cluster coherence k=10:", coherence_k10, "\n")

coh_df <- data.frame(k = k_range, coherence = coherence_scores)
ggplot(coh_df, aes(x = k, y = coherence)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = 10, col = "red", linetype = "dashed") +
  labs(x = "Number of clusters (k)", y = "Mean within-cluster Jaccard similarity",
       title = "Within-cluster coherence across k") +
  theme_classic()
ggsave("figures/coherence_scores.pdf", width = 6, height = 4)

# Summary
cat("ARI (spectral vs hierarchical):", ari, "\n")
cat("Mean silhouette:", mean(sil_k10[, 3]), "\n")
cat("Within-cluster Jaccard coherence:", coherence_k10, "\n")