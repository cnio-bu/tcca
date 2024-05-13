library(beyondcell)
library(tidyverse)
library(ComplexHeatmap)

## function fix for factoextra
my_fviz_nbclust <- function(x, print.summary = TRUE, barfill = "steelblue", barcolor = "steelblue"){
    best_nc <- x$Best.nc
    best_nc <- as.data.frame(t(best_nc), stringsAsFactors = TRUE)
    best_nc$Number_clusters <- as.factor(best_nc$Number_clusters)
    
    ss <- summary(best_nc$Number_clusters)
    cat("Among all indices: \n===================\n")
    for (i in 1:length(ss)) {
        cat("*", ss[i], "proposed ", names(ss)[i], "as the best number of clusters\n")
    }
    cat("\nConclusion\n=========================\n")
    cat("* According to the majority rule, the best number of clusters is ", 
        names(which.max(ss)), ".\n\n")
    
    df <- data.frame(Number_clusters = names(ss), freq = ss, 
                     stringsAsFactors = TRUE)
    p <- ggpubr::ggbarplot(df, x = "Number_clusters", y = "freq", 
                           fill = "steelblue", color = "steelblue") +
        ggplot2::labs(x = "Number of clusters k", 
                      y = "Frequency among all indices",
                      title = paste0("Optimal number of clusters - k = ", 
                                     names(which.max(ss))))
    p
}

bcs <- readRDS("results/single_sample/breast_sunny_wu.rds")
single_sample_bc <- bcs[[6]]

mat <- single_sample_bc@normalized

rownames(mat)

## make a zscore of the full mat to ease the PCA. Each drug value is scaled over cells
mat2 <- scale(x = t(mat), center = TRUE, scale = TRUE)
mat2 <- t(mat2)

## first, let's see value distributions of the most variable features
top_rv <- matrixStats::rowVars(mat2)
top_rv <- sort(top_rv, decreasing = TRUE)

top10_drugs <- head(top_rv, n = 10)

top10_dt <- as.data.frame(t(mat2[names(top10_drugs), ]))
top10_dt_long <- top10_dt %>%
    pivot_longer(names_to = "drug", values_to = "enrichment", cols = everything())

dist_top_drugs <- ggplot(data = top10_dt_long, aes(x = enrichment)) +
    scale_x_continuous(name = "zscore BCS", n.breaks = 10) +
    scale_y_continuous(name = "Density", n.breaks = 10) +
    geom_rug(alpha = 0.2) +
    geom_histogram(bins = 50) +
    ggtitle("BCS distribution of most variable drugs in a single breast cancer patient") +
    theme_bw() +
    theme(
        text = element_text(family = "Arial")
    ) +
    facet_wrap(~drug, nrow = 5, ncol = 2)


ggsave(
    filename = "results/single_sample/clustering_expr/top_10_variable_drugs_dist.png",
    plot = dist_top_drugs,
    dpi = 300,
    height = 7,
    width = 14
)

## Perform PCA
missing_drugs <- rowSums(mat2)
missing_drugs <- missing_drugs[is.na(missing_drugs)]

mat2 <- mat2[!(rownames(mat2) %in% names(missing_drugs)), ]

pra <- prcomp(x = t(mat2), center = FALSE, scale. = FALSE)
rotation_pca <- pra$x

pca_plot <- ggplot(data = as.data.frame(rotation_pca), aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.2) +
    scale_x_continuous(name = "PC1", n.breaks = 10) +
    scale_y_continuous(name = "PC2", n.breaks = 10) +
    ggtitle("PCA of the scaled BCS for 70 drugs and 1237 cells") +
    theme_bw() +
    theme(text = element_text(family =  "Arial"))

pca_plot_23 <- ggplot(data = as.data.frame(rotation_pca), aes( x = PC2, y = PC3)) +
    geom_point(alpha = 0.2) +
    scale_x_continuous(name = "PC2", n.breaks = 10) +
    scale_y_continuous(name = "PC3", n.breaks = 10) +
    ggtitle("PCA of the scaled BCS for 70 drugs and 1237 cells") +
    theme_bw() +
    theme(text = element_text(family =  "Arial"))


ggsave(
    filename = "results/single_sample/clustering_expr/pca_of_90drugs.png",
    plot = pca_plot,
    dpi = 100,
    height = 7,
    width = 7
)

ggsave(
    filename = "results/single_sample/clustering_expr/pca_of_90drugs_2_3.png",
    plot = pca_plot_23,
    dpi = 100,
    height = 7,
    width = 7
)

# the contribution to the total variance for each component
percentVar <- data.frame("explained_var" = (pra$sdev^2 / sum(pra$sdev^2 )))
percentVar$component <- paste0("PC", c(1:70))

percentVar$component <- as.factor(percentVar$component)
percentVar$component <- fct_reorder(percentVar$component, percentVar$explained_var, .desc = TRUE)

explained_var <- ggplot(data = percentVar[1:10, ], aes(x = component, y = explained_var)) +
    geom_bar(stat = "identity") + 
    scale_x_discrete(name = "") +
    scale_y_continuous("% total variance explained", labels = scales::percent,  n.breaks = 10) +
    theme_bw() +
    theme(text = element_text(family = "Arial"))

ggsave(
    filename = "results/single_sample/clustering_expr/explained_variance_10_components.png",
    plot = explained_var,
    dpi = 100,
    height = 7,
    width = 7
)

mat3 <- t(apply(mat, 1, scales::rescale, to = c(0, 1)))
mat3 <- round(mat3, digits = 2)

top_rv <- matrixStats::rowVars(mat3)
top_rv <- sort(top_rv, decreasing = TRUE)

top10_drugs <- head(top_rv, n = 10)

top10_dt <- as.data.frame(t(mat3[names(top10_drugs), ]))
top10_dt_long <- top10_dt %>%
    pivot_longer(names_to = "drug", values_to = "enrichment", cols = everything())


dist_top_drugs_scaled_2 <- ggplot(data = top10_dt_long, aes(x = enrichment)) +
    scale_x_continuous(name = "BCS", n.breaks = 10) +
    scale_y_continuous(name = "Density", n.breaks = 10) +
    geom_rug(alpha = 0.2) +
    geom_histogram(bins = 50) +
    ggtitle("BCS distribution of most variable drugs in a single breast cancer patient with BC scaler") +
    theme_bw() +
    theme(
        text = element_text(family = "Arial")
    ) +
    facet_wrap(~drug, nrow = 5, ncol = 2)


ggsave(
    filename = "results/single_sample/clustering_expr/top10_variable_drugs_bcscaler.png",
    plot = dist_top_drugs_scaled_2,
    dpi = 300,
    height = 7,
    width = 7
)

## redo PCA
pra <- prcomp(x = t(mat3), center = FALSE, scale. = FALSE)
rotation_pca <- pra$x

pca_plot <- ggplot(data = as.data.frame(rotation_pca), aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.2) +
    scale_x_continuous(name = "PC1", n.breaks = 10) +
    scale_y_continuous(name = "PC2", n.breaks = 10) +
    ggtitle("PCA of the scaled BCS for 70 drugs and 1237 cells using BC scaler") +
    theme_bw() +
    theme(text = element_text(family =  "Arial"))

pca_plot_23 <- ggplot(data = as.data.frame(rotation_pca), aes(x = PC2, y = PC3)) +
    geom_point(alpha = 0.2) +
    scale_x_continuous(name = "PC2", n.breaks = 10) +
    scale_y_continuous(name = "PC3", n.breaks = 10) +
    ggtitle("PCA of the scaled BCS for 70 drugs and 1237 cellss using BC scaler") +
    theme_bw() +
    theme(text = element_text(family =  "Arial"))

ggsave(
    filename = "results/single_sample/clustering_expr/pca_of_90drugs_scaledbc.png",
    plot = pca_plot,
    dpi = 100,
    height = 7,
    width = 7
)

ggsave(
    filename = "results/single_sample/clustering_expr/pca_of_90drugs_2_3_scaledbc.png",
    plot = pca_plot_23,
    dpi = 100,
    height = 7,
    width = 7
)

## HEATMAP of the zscaled bcscores
png("results/single_sample/clustering_expr/heatmap_scaled_sig.png", width = 9, height = 12, res = 100, units = "in")
heat <- ComplexHeatmap::Heatmap(
    matrix = mat2,
    name = "Zscaled bcscores",
    #col = c("blue", "white", "red"),
    column_split = 10,
    row_split = 10,
    cluster_row_slices = TRUE,
    cluster_column_slices = TRUE,
    row_names_gp = gpar(fontsize = 6),
    show_column_names = FALSE
)
draw(heat)
dev.off()

cor.mat <- cor(t(mat2))
png("results/single_sample/clustering_expr/heatmap_scaled_sig_correlations.png", width = 9, height = 9, res = 100, units = "in")
heat2 <- ComplexHeatmap::Heatmap(
    matrix = cor.mat,
    name = "Pearson correlation",
    #col = c("blue", "white", "red"),
    column_split = 4,
    row_split = 4,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    cluster_row_slices = TRUE,
    cluster_column_slices = TRUE,
    row_names_gp = gpar(fontsize = 6),
    show_column_names = FALSE
)
draw(heat2)
dev.off()

cor.mat2 <- cor(mat2)
png("results/single_sample/clustering_expr/heatmap_scaled_cell_correlations.png", width = 9, height = 9, res = 100, units = "in")
heat3 <- ComplexHeatmap::Heatmap(
    matrix = cor.mat2,
    show_row_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    cluster_row_slices = TRUE,
    cluster_column_slices = TRUE,
    show_column_names = FALSE
)
draw(heat3)
dev.off()


library(Nbclust)
library(factoextra)

test <- NbClust::NbClust(
    data = t(mat2),
    method = "kmeans",
    min.nc = 2,
    max.nc = 10
    )

a <- my_fviz_nbclust(test)

test2 <- NbClust::NbClust(
    data = t(mat3),
    method = "kmeans",
    min.nc = 2,
    max.nc = 10
)

ggsave(plot = a, filename = "results/single_sample/clustering_expr/best_clusters_bc_zscaled.png")


single_sample_bc <- bcUMAP(bc = single_sample_bc, pc = 10, k.neighbors = 20, res = 0.2)
bcclusters <- bcClusters(single_sample_bc, idents = "bc_clusters_res.0.2")

ggsave(plot = bcclusters, filename = "results/single_sample/clustering_expr/therapeutic_clusters_umap_with_bcscaled.png")

single_sample_bc@scaled <- mat2
single_sample_bc <- bcUMAP(bc = single_sample_bc, pc = 10, k.neighbors = 20, res = 0.2)


raw_uwot <- uwot::umap(X  = prcomp(t(mat2))$x, seed = 42)
raw_uwot <- as.data.frame(raw_uwot)
colnames(raw_uwot) <- c("UMAP 1", "UMAP 2")



raw_uwot$cell <- rownames(raw_uwot)

bc_clusts <- single_sample_bc@meta.data
bc_clusts$cell <- rownames(bc_clusts)

raw_uwot_annot <- raw_uwot %>%
    left_join(bc_clusts[, c("cell", "bc_clusters_res.0.2")], by = "cell")

uwot_scaled <- ggplot(data = raw_uwot_annot, aes(x = `UMAP 1`, y = `UMAP 2`, color = bc_clusters_res.0.2)) +
    geom_point() +
    theme_minimal()

ggsave(plot = uwot_scaled, filename = "results/single_sample/uwot_bcscores_zscaled.png", dpi = 300)


pra_annot <- pra$x
pra_annot <- as.data.frame(pra_annot)
pra_annot$cell <- rownames(pra_annot)

pra_annot <- pra_annot[, c("cell", "PC1", "PC2", "PC3")]

pra_annot <- pra_annot %>%
    left_join(bc_clusts[, c("cell", "bc_clusters_res.0.2")], by = "cell")

pra_annot_scaled <- ggplot(data = pra_annot, aes(x = PC1, y = PC2, color = bc_clusters_res.0.2)) +
    geom_point() +
    theme_minimal()

pra_annot_scaled_23 <- ggplot(data = pra_annot, aes(x = PC2, y = PC3, color = bc_clusters_res.0.2)) +
    geom_point() +
    theme_minimal()

ggsave(plot = pra_annot_scaled, filename = "results/single_sample/clustering_expr/pca_of90drugs_scaledbc_with_tcs.png")
ggsave(plot = pra_annot_scaled_23, filename = "results/single_sample/clustering_expr/pca_of90drugs_scaledbc_with_tcs_23.png")


pra <- prcomp(x = t(mat2), center = FALSE, scale. = FALSE)
rotation_pca <- pra$x


rotation_pca <- as.data.frame(rotation_pca)
rotation_pca$cell <- rownames(rotation_pca)

rotation_pca <- rotation_pca[, c("cell", "PC1", "PC2", "PC3")]

rotation_pca <- rotation_pca %>%
    left_join(bc_clusts[, c("cell", "bc_clusters_res.0.2")], by = "cell")

pra_annot_scaled <- ggplot(data = rotation_pca, aes(x = PC1, y = PC2, color = bc_clusters_res.0.2)) +
    geom_point() +
    theme_minimal()

pra_annot_scaled_23 <- ggplot(data = rotation_pca, aes(x = PC2, y = PC3, color = bc_clusters_res.0.2)) +
    geom_point() +
    theme_minimal()


ggsave(plot = pra_annot_scaled, filename = "results/single_sample/clustering_expr/pca_of90drugs_zscore_with_tcs.png")
ggsave(plot = pra_annot_scaled_23, filename = "results/single_sample/clustering_expr/pca_of90drugs_zscore_with_tcs_23.png")


