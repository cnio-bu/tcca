library(dplyr)
library(tidyverse)
library(cluster)
library(mclust)
library(kernlab)
library(fpc)
library(parallel)

set.seed(123)

# Paths
base_dir <- "/storage/scratch01/shared/projects/bc-meta/single_cell/sctherapy"
out_dir <- file.path(base_dir, "results/clustering_metrics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")

# Load data
cat("[1/6] Loading data...\n")
similarity_matrix <- readRDS(file.path(base_dir, "results/jaccard_matrix.rds"))
cluster_assignment <- readRDS(file.path(base_dir, "results/speclustering_reordered.rds"))
dist_matrix <- as.dist(1 - similarity_matrix)
k_range <- 2:14
n_cores <- max(1, detectCores() - 1)

cat("Matrix dimensions:", dim(similarity_matrix), "\n")
cat("N subclones:", length(cluster_assignment), "\n")
cat("N cores available:", n_cores, "\n")

#  1. ARI spectral vs hierarchical
cat("[2/6] ARI spectral vs hierarchical...\n")
hc <- hclust(dist_matrix, method = "ward.D2")
hc_clusters <- cutree(hc, k = 10)
ari <- adjustedRandIndex(as.integer(cluster_assignment), hc_clusters)
cat("ARI spectral vs hierarchical:", round(ari, 3), "\n")

# 2. Bootstrap ARI (parallelized)
cat("[3/6] Bootstrap ARI (100 runs, parallelized)...\n")
ari_runs <- unlist(mclapply(1:100, function(i) {
    set.seed(i)
    run_i <- specc(similarity_matrix, centers = 10)
    adjustedRandIndex(as.integer(cluster_assignment), as.integer(run_i))
}, mc.cores = n_cores))

cat(
    "Bootstrap ARI mean ± sd:",
    round(mean(ari_runs), 3), "±", round(sd(ari_runs), 3), "\n"
)

saveRDS(ari_runs, file.path(out_dir, "ari_bootstrap_runs.rds"))

# 3. Clusterboot Jaccard stability
cat("[4/6] Clusterboot Jaccard stability...\n")

specc_CBI <- function(x, k, ...) {
    result <- specc(x, centers = k)
    clusters <- as.integer(result)
    list(
        result = result,
        nc = k,
        clusterlist = lapply(1:k, function(cl) clusters == cl),
        partition = clusters,
        clustermethod = "specc"
    )
}

stab <- clusterboot(
    similarity_matrix,
    B = 100,
    clustermethod = specc_CBI,
    k = 10,
    seed = 42,
    count = TRUE
)

saveRDS(stab, file.path(out_dir, "clusterboot_stability.rds"))

# Stability table per cluster
stab_df <- data.frame(
    TC = paste0("TC", 1:10),
    jaccard = round(stab$bootmean, 3),
    dissolved = stab$bootbrd,
    recovered = stab$bootrecover,
    stability = case_when(
        stab$bootmean >= 0.75 ~ "stable",
        stab$bootmean >= 0.60 ~ "borderline",
        TRUE ~ "unstable"
    )
) %>% arrange(jaccard)

write.table(stab_df,
    file.path(out_dir, "clusterboot_per_TC.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)
print(stab_df)

# 4. Silhouette and coherence across k
cat("[5/6] Silhouette and coherence across k...\n")

compute_coherence <- function(clusters, sim_matrix) {
    mean(sapply(unique(clusters), function(cl) {
        idx <- which(clusters == cl)
        if (length(idx) < 2) {
            return(NA)
        }
        vals <- sim_matrix[idx, idx]
        mean(vals[lower.tri(vals)])
    }), na.rm = TRUE)
}

# Run spectral clustering for silhoutte and coherence score computation
specc_results <- lapply(k_range, function(k) {
    cat("  Running specc k =", k, "\n")
    set.seed(123)
    specc(similarity_matrix, centers = k)
})
names(specc_results) <- paste0("k", k_range)

silhouette_scores <- sapply(specc_results, function(res) {
    sil <- silhouette(as.integer(res), dist_matrix)
    mean(sil[, 3])
})

coherence_scores <- sapply(specc_results, function(res) {
    compute_coherence(as.integer(res), similarity_matrix)
})

# k=10 with clustering saved (specc is not recalculated)
sil_k10 <- silhouette(as.integer(cluster_assignment), dist_matrix)
coherence_k10 <- compute_coherence(as.integer(cluster_assignment), similarity_matrix)
cat("Silhouette k=10:", round(mean(sil_k10[, 3]), 3), "\n")
cat("Within-cluster coherence k=10:", round(coherence_k10, 3), "\n")

# 5. Plots
cat("[6/6] Generating plots...\n")

# Plot A: silhouette + coherence across k
metrics_df <- data.frame(
    k = k_range,
    Silhouette = silhouette_scores,
    Coherence = coherence_scores
) %>%
    pivot_longer(-k, names_to = "metric", values_to = "value")

p_metrics <- ggplot(metrics_df, aes(x = k, y = value, color = metric)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_vline(xintercept = 10, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_color_manual(
        values = c("Silhouette" = "#4575b4", "Coherence" = "#d73027"),
        labels = c("Silhouette score", "Within-cluster Jaccard")
    ) +
    scale_x_continuous(breaks = k_range) +
    facet_wrap(~metric,
        scales = "free_y", ncol = 1,
        labeller = labeller(metric = c(
            "Silhouette" = "Mean silhouette score",
            "Coherence" = "Mean within-cluster Jaccard similarity"
        ))
    ) +
    labs(
        x = "Number of clusters (k)", y = NULL,
        title = "Cluster quality metrics across k",
        subtitle = paste0(
            "ARI spectral vs hierarchical = ", round(ari, 3),
            "  |  Bootstrap ARI = ", round(mean(ari_runs), 3),
            " ± ", round(sd(ari_runs), 3)
        )
    ) +
    theme_classic(base_size = 11) +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "gray40")
    )

ggsave(file.path(out_dir, "cluster_metrics_across_k.pdf"),
    p_metrics,
    width = 5, height = 6
)

# Plot B — per-cluster silhouette
sil_per_cluster <- as.data.frame(sil_k10[, ]) %>%
    mutate(TC = paste0("TC", cluster)) %>%
    group_by(TC) %>%
    summarise(mean_sil = mean(sil_width), n = n()) %>%
    arrange(mean_sil)

p_sil_cluster <- ggplot(
    sil_per_cluster,
    aes(x = reorder(TC, mean_sil), y = mean_sil, fill = mean_sil)
) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_text(aes(label = paste0("n=", n)),
        hjust = -0.2, size = 3, color = "gray30"
    ) +
    scale_fill_gradient2(
        low = "#d73027", mid = "gray90", high = "#4575b4",
        midpoint = 0, guide = "none"
    ) +
    coord_flip() +
    expand_limits(y = max(sil_per_cluster$mean_sil) * 1.15) +
    labs(
        x = NULL, y = "Mean silhouette score",
        title = "Per-cluster silhouette score (k=10)",
        subtitle = "Clusters with score < 0 are candidates for merging/removal"
    ) +
    theme_classic(base_size = 11)

ggsave(file.path(out_dir, "silhouette_per_cluster.pdf"),
    p_sil_cluster,
    width = 5, height = 4
)

# Plot C - clusterboot Jaccard per TC
p_stab <- ggplot(
    stab_df,
    aes(x = reorder(TC, jaccard), y = jaccard, fill = stability)
) +
    geom_col() +
    geom_hline(yintercept = 0.75, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0.60, linetype = "dotted", color = "gray50") +
    scale_fill_manual(values = c(
        "stable" = "#4575b4",
        "borderline" = "#fee090",
        "unstable" = "#d73027"
    )) +
    coord_flip() +
    labs(
        x = NULL, y = "Mean Jaccard similarity (bootstrap)",
        title = "Cluster stability (clusterboot, B=100)",
        subtitle = "Dashed = 0.75 (stable), dotted = 0.60 (borderline)"
    ) +
    theme_classic(base_size = 11)

ggsave(file.path(out_dir, "clusterboot_stability.pdf"),
    p_stab,
    width = 5, height = 4
)

# 6. Summary table
summary_df <- data.frame(
    Metric = c(
        "ARI spectral vs hierarchical",
        "Bootstrap ARI mean",
        "Bootstrap ARI sd",
        "Mean silhouette (k=10)",
        "Within-cluster Jaccard coherence (k=10)"
    ),
    Value = round(c(
        ari,
        mean(ari_runs),
        sd(ari_runs),
        mean(sil_k10[, 3]),
        coherence_k10
    ), 3)
)

write.table(summary_df,
    file.path(out_dir, "cluster_stability_summary.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

cat("\n── Summary ───────────────────────────────────────────\n")
print(summary_df)
cat("\nAll outputs saved to:", out_dir, "\n")