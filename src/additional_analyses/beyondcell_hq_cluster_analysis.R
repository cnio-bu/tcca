library(BPCells)
library(matrixStats)
library(tidyverse)

## We perform clustering exploratory analysis
## using a sketch of 50k cells x 90 signatures

mat <- BPCells::open_matrix_dir(dir = "results/beyondcell_bp/sketch_mat_beyondcell")
mat <- as.matrix(mat)


## first, let's see value distributions of the most variable features
top_rv <- matrixStats::rowVars(mat)
top_rv <- sort(top_rv, decreasing = TRUE)

top10_drugs <- head(top_rv, n = 10)

top10_dt <- as.data.frame(t(mat[names(top10_drugs), ]))
top10_dt_long <- top10_dt %>%
    pivot_longer(names_to = "drug", values_to = "enrichment", cols = everything())


dist_top_drugs <- ggplot(data = top10_dt_long, aes(x = enrichment)) +
    scale_x_continuous(name = "BCS", n.breaks = 10) +
    scale_y_continuous(name = "Density", n.breaks = 10) +
    geom_rug(alpha = 0.2) +
    geom_histogram(bins = 50) +
    ggtitle("BCS distribution of most variable drugs in the 50k. sketch") +
    theme_bw() +
    theme(
        text = element_text(family = "Arial")
    ) +
    facet_wrap(~drug, nrow = 5, ncol = 2)

ggsave(
    filename = "results/tcca/clustering_expr/top_10_variable_drugs_dist.png",
    plot = dist_top_drugs,
    dpi = 300,
    height = 7,
    width = 14
    )

## make a zscore of the full mat to ease the PCA. Each drug value is scaled over cells
mat <- scale(x = t(mat), center = TRUE, scale = TRUE)
mat <- t(mat)

top10_dt_scaled <- as.data.frame(t(mat[names(top10_drugs), ]))
top10_dt_long_scaled <- top10_dt_scaled %>%
    pivot_longer(names_to = "drug", values_to = "enrichment", cols = everything())

dist_top_drugs_scaled <- ggplot(data = top10_dt_long_scaled, aes(x = enrichment)) +
    scale_x_continuous(name = "Scaled BCS", n.breaks = 10) +
    scale_y_continuous(name = "Density", n.breaks = 10) +
    geom_rug(alpha = 0.2) +
    geom_histogram(bins = 50) +
    ggtitle("BCS distribution of most variable drugs in the 50k. sketch") +
    theme_bw() +
    theme(
        text = element_text(family = "Arial")
    ) +
    facet_wrap(~drug, nrow = 5, ncol = 2)

ggsave(
    filename = "results/tcca/clustering_expr/top_10_variable_drugs_dist_scaled.png",
    plot = dist_top_drugs_scaled,
    dpi = 300,
    height = 7,
    width = 14
)

## Perform PCA
pra <- prcomp(x = t(mat), center = FALSE, scale. = FALSE)
rotation_pca <- pra$x

pca_plot <- ggplot(data = as.data.frame(rotation_pca), aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.2) +
    scale_x_continuous(name = "PC1", n.breaks = 10) +
    scale_y_continuous(name = "PC2", n.breaks = 10) +
    ggtitle("PCA of the scaled BCS for 90 drugs and 50.000 cells") +
    theme_bw() +
    theme(text = element_text(family =  "Arial"))

pca_plot_23 <- ggplot(data = as.data.frame(rotation_pca), aes( x = PC2, y = PC3)) +
    geom_point(alpha = 0.2) +
    scale_x_continuous(name = "PC2", n.breaks = 10) +
    scale_y_continuous(name = "PC3", n.breaks = 10) +
    ggtitle("PCA of the scaled BCS for 90 drugs and 50.000 cells") +
    theme_bw() +
    theme(text = element_text(family =  "Arial"))


ggsave(
    filename = "results/tcca/clustering_expr/pca_of_90drugs.png",
    plot = pca_plot,
    dpi = 100,
    height = 7,
    width = 7
    )

ggsave(
    filename = "results/tcca/clustering_expr/pca_of_90drugs_2_3.png",
    plot = pca_plot_23,
    dpi = 100,
    height = 7,
    width = 7
)

# the contribution to the total variance for each component
percentVar <- data.frame("explained_var" = (pra$sdev^2 / sum(pra$sdev^2 )))
percentVar$component <- paste0("PC", c(1:90))

percentVar$component <- as.factor(percentVar$component)
percentVar$component <- fct_reorder(percentVar$component, percentVar$explained_var, .desc = TRUE)

explained_var <- ggplot(data = percentVar[1:10, ], aes(x = component, y = explained_var)) +
    geom_bar(stat = "identity") + 
    scale_x_discrete(name = "") +
    scale_y_continuous("% total variance explained", labels = scales::percent,  n.breaks = 10) +
    theme_bw() +
    theme(text = element_text(family = "Arial"))

ggsave(
    filename = "results/tcca/clustering_expr/explained_variance_10_components.png",
    plot = explained_var,
    dpi = 100,
    height = 7,
    width = 7
)


## Partition cluster optimization
library(dbscan)
library(NbClust)

hdb <- hdbscan(
    x = t(mat),
    minPts = 100,
    verbose = TRUE,
    gen_hdbscan_tree = TRUE,
    gen_simplified_tree = TRUE
    )


png(
    filename = "results/tcca/clustering_expr/simplified_hdb_scan.png",
    width = 7,
    height = 7,
    units = "in",
    res = 100
    )

plot(hdb)
dev.off()

plot(hdb$hc, main="HDBSCAN* Hierarchy")
