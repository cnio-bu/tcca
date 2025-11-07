library(Seurat)
library(BPCells)
library(UCell)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

setwd("/storage/scratch01/shared/projects/bc-meta/functional_nmf")

# # Load seurat object with malignant cells
# seu_lvl2 <- readRDS("../single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
# malignant <- subset(seu_lvl2, subset = malignancy == TRUE)

# # Load the list of metaprograms
# mp_list <- readRDS("sample_wise/metaprograms_cpm/mp_list_reordered.rds")
# names(mp_list) <- paste0("MP", 1:length(mp_list))
# mp_ucell <- AddModuleScore_UCell(malignant, features = mp_list)

# # Save seurat object with UCell scores in metadata
# saveRDS(mp_ucell, "sample_wise/seurat_mps_ucell.rds")
# # Save table of UCell scores
# write.table(mp_ucell@meta.data[, paste0("MP", 1: 43, "_UCell")], 
#             "sample_wise/mps_ucell_scores.tsv")

######################## PLOT UCELL SCORES PER METAPROGRAM #####################
seu_mp <- readRDS("sample_wise/seurat_mps_ucell.rds")

# Add clonal information and sctherapy clusters to the metadata
seu_subclones <- readRDS("../single_cell/sctherapy/results/seu_subclones.rds")

seu_mp <- subset(seu_mp, cells = colnames(seu_subclones))
new_columns <- setdiff(colnames(seu_subclones@meta.data), colnames(seu_mp@meta.data))
new_metadata <- seu_subclones@meta.data[new_columns] %>%
    rownames_to_column(var = "cell_id")
seu_mp@meta.data <- seu_mp@meta.data %>%
    rownames_to_column(var = "cell_id") %>%
    left_join(new_metadata, by = "cell_id") %>%
    column_to_rownames(var = "cell_id")

# Create a sketch of 5000 cells
seu_mp <- NormalizeData(seu_mp, 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
    )

seu_mp <- FindVariableFeatures(seu_mp, selection.method = "vst", nfeatures = 2000)
seu_mp <-  ScaleData(seu_mp, features = rownames(seu_mp))
hvg <- VariableFeatures(seu_mp)
seu_mp <- JoinLayers(seu_mp)
seu_mp <- SketchData(
  object = seu_mp,
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch",
  features = hvg
)

saveRDS(seu_mp, "sample_wise/seu_mps_sketch.rds")


### 1. Plot a heatmap of raw UCell scores for a sketch of 5k cells
# Subset 5k cells with their ucell scores
cells <- colnames(seu_mp[["sketch"]]$counts)
ucell_mat_subset <- as.matrix(t(seu_mp@meta.data[
    cells,
    grepl("UCell", colnames(seu_mp@meta.data))
]))

# Create top heatmap annotations (sample site and cluster)
source("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")
translat_human_sites <- c(
    "adrenal_gland" = "Adrenal gland",
    "bladder" = "Bladder",
    "bone_marrow" = "Bone marrow",
    "brain" = "Brain",
    "breast" = "Breast",
    "colon" = "Colon",
    "skin" = "Skin",
    "esophagus" = "Esophagus",
    "oesophagus" = "Esophagus",
    "kidney" = "Kidney",
    "liver" = "Liver",
    "lung" = "Lung",
    "lymph_node" = "Lymph node",
    "other" = "Other",
    "ovary" = "Ovary",
    "pancreas" = "Pancreas",
    "prostate" = "Prostate",
    "soft_tissue" = "Soft tissue"
)

top_annotation_df <- seu_mp@meta.data %>%
    filter(cell %in% cells) %>%
    select(tumor_site, clusters) %>%
    mutate(tumor_site = case_when(
        tumor_site %in% names(translat_human_sites) ~ tumor_site,
        TRUE ~ "other"
    )) %>%
    arrange(clusters, tumor_site)

top_annotation_df$tumor_site <- translat_human_sites[
    top_annotation_df$tumor_site
]

colnames(top_annotation_df) <- c(
    "Sample site",
    "scTherapy cluster"
)

pals <- list(
    "Sample site" = tumor_sites_colors,
    "scTherapy cluster" = sctherapy_colors
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = top_annotation_df,
    which = "column",
    col = pals,
    annotation_name_side = "left",
    annotation_name_rot = 0,
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
    annotation_legend_param = list(
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 12),
        title_gap = unit(10, "mm")
    ),
    show_legend = c(
        "Sample site" = FALSE, "scTherapy cluster" = FALSE
    )
)

# Customize legends for tumor site and cluster
tumor_site_legend <- Legend(
    at = names(pals$`Sample site`),
    legend_gp = gpar(fill = pals$`Sample site`),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    ncol = 3, # Split Group legend into 2 columns
    gap = unit(5, "mm"),
    title = "Sample site"
)

cluster_legend <- Legend(
    at = names(pals$`scTherapy cluster`),
    legend_gp = gpar(fill = pals$`scTherapy cluster`),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    ncol = 2, # Split Group legend into 2 columns
    gap = unit(5, "mm"),
    title = "scTherapy cluster"
)


col_fun <- colorRamp2(
    breaks = c(0, 0.125, 0.25, 0.375, 0.5), # set appropriate min/max
    colors = c("#3B4CC0", "#78D0AA", "#F7F7BD", "#F89560", "#B8122A")
)

# Create heatmap
heat <- ComplexHeatmap::Heatmap(
    mat = ucell_mat_subset[, rownames(top_annotation_df)],
    col = col_fun,
    top_annotation = top_annotation,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = TRUE,
    show_column_dend = FALSE,
    column_split = top_annotation_df$`scTherapy cluster`,
    clustering_distance_columns = "pearson",
    clustering_distance_rows = "pearson",
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 12),
    row_names_side = "right",
    row_title = "Metaprograms",
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_title = "Sketch of 5000 malignant cells",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_title_side = "left",
    column_title_side = "bottom",
    heatmap_legend_param = list(
        at = seq(0, 0.5, by = 0.1),
        title = "UCell score",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 12),
        title_gap = unit(10, "mm")
    ),
    heatmap_width = unit(10, "in"),
    heatmap_height = unit(8, "in"),
    use_raster = TRUE,
    raster_quality = 5
)

png(
    file = "sample_wise/figures/heatmap_ucell.png",
    res = 500,
    width = 13,
    height = 12,
    units = "in"
)
draw(heat,
    annotation_legend_side = "top",
    heatmap_legend_side = "right",
    annotation_legend_list = list(tumor_site_legend, cluster_legend)
)
dev.off()




### 2. Dot plot of mean UCell scores per Metaprogram
# Subset UCell scores per metaprogram
ucell_df <- seu_mp@meta.data %>%
    select(clusters, contains("UCell")) %>%
    select(-CIN70_UCell) %>%
    pivot_longer(-clusters, names_to = "metaprogram", values_to = "score")

# Compute per-signature 75th percentile
mp_thresholds <- ucell_df %>%
    group_by(metaprogram) %>%
    summarise(thresh_75 = quantile(score, 0.75, na.rm = TRUE))

# Join thresholds and compute mean + % above threshold
summary_df <- ucell_df %>%
    left_join(mp_thresholds, by = "metaprogram") %>%
    group_by(clusters, metaprogram) %>%
    summarise(
        mean_score = mean(score, na.rm = TRUE),
        pct_active = mean(score > thresh_75, na.rm = TRUE) * 100,
        .groups = "drop"
    ) %>%
    mutate(
        metaprogram = gsub("_UCell", "", metaprogram),
        metaprogram = factor(metaprogram, levels = paste0("MP", 1:43)),
        clusters = factor(clusters, levels = as.character(10:1))
        )

bubble_mean <- ggplot(summary_df, aes(x = metaprogram, y =  clusters)) +
    geom_point(aes(size = pct_active, color = mean_score)) +
    scale_color_gradientn(
        colors = c("#3B4CC0", "#78D0AA", "#F7F7BD", "#F89560", "#B8122A"),
        limits = c(0, 0.5),
        labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
        oob = scales::squish
    ) +
    scale_size(range = c(1, 6)) +
    theme_minimal() +
    labs(
        color = "Mean UCell score",
        size = "% High-scoring cells\n(>75th pct)",
        x = "Metaprogram", y = "Cluster"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(face = "bold")
    )

ggsave("sample_wise/figures/bubble_mean05.png",
    plot = bubble_mean, 
    width = 16, 
    height = 5
)

### 3. Dot plot of scaled UCell scores per Metaprogram
# z-scaled ucell scores
ucell_df_scaled <- ucell_df %>%
    group_by(metaprogram) %>%
    mutate(zscore = as.numeric(scale(score))) %>%
    ungroup()

# Compute mean of z-scaled ucell scores
summary_df <- ucell_df_scaled %>%
    left_join(mp_thresholds, by = "metaprogram") %>%
    group_by(clusters, metaprogram) %>%
    summarise(
        mean_zscore = mean(zscore, na.rm = TRUE),
        pct_active = mean(score > thresh_75, na.rm = TRUE) * 100,
        .groups = "drop"
    ) %>%
    mutate(
        metaprogram = gsub("_UCell", "", metaprogram),
        metaprogram = factor(metaprogram, levels = paste0("MP", 1:43)),
        clusters = factor(clusters, levels = as.character(10:1))
        )

# Scale mean of z-scaled ucell score between 0 and 1 for clarity
summary_df <- summary_df %>%
    mutate(
        zscore_rescaled = (mean_zscore - min(mean_zscore, na.rm = TRUE)) /
            (max(mean_zscore, na.rm = TRUE) - min(mean_zscore, na.rm = TRUE))
    )

bubble_mean <- ggplot(summary_df, aes(x = metaprogram, y = clusters)) +
    geom_point(aes(size = pct_active, color = mean_zscore)) +
    scale_color_gradientn(
        colors = c("#3B4CC0", "#78D0AA", "#F7F7BD", "#F89560", "#B8122A"),
        limits = c(-0.5, 1),
        breaks = seq(-0.5, 1, by = 0.3),
        oob = scales::squish
    ) +
    scale_size(range = c(1, 6)) +
    theme_minimal() +
    labs(
        color = "Mean UCell\n(z-score)",
        size = "% High-scoring cells\n(>75th pct)",
        x = "Metaprogram", y = "Cluster"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12)
    )

ggsave("sample_wise/figures/bubble_mean_zscaled.png",
    plot = bubble_mean,
    width = 15,
    height = 4
)

pdf("sample_wise/figures/bubble_mean_zscaled.pdf",
    width = 15,
    height = 4
)
bubble_mean
dev.off()



### Compute CIN70 signature enrichment per TC
cin70_sig <- read.xlsx("sample_wise/cin70.xlsx",
    sheet = "Table S3",
    startRow = 2
)

cin70_sig <- cin70_sig$CIN70.signature
seu_mp <- AddModuleScore_UCell(seu_mp, features = list(CIN70 = cin70_sig))


vl <- ggplot(seu_mp@meta.data, aes(x = clusters, y = CIN70_UCell)) +
    geom_violin(aes(fill = clusters), scale = "width", trim = FALSE) +
    scale_fill_manual(values = sctherapy_colors) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    theme_bw(base_size = 9) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12)
    ) +
    labs(x = "Cluster", y = "CIN70 score") +
    NoLegend()

ggsave("sample_wise/figures/vln_cin70.png", plot = vl, width = 5, height = 5, dpi = 300)

pdf("sample_wise/figures/vln_cin70.pdf",
    width = 5,
    height = 4
)
vl
dev.off()


### Compute transcriptomic heterogeneity per TC
library(entropy)
# Get log-normalized expression matrix
expr_mat <- as.matrix(seu_mp[["RNA"]]$data)

# Discretize (bin) expression to stabilize entropy calculation
# You can use 10 bins (can be tuned)
entropy_per_cell <- apply(expr_mat, 2, function(x) {
    entropy(discretize(x, numBins = 10))
})

# Add to Seurat metadata
seu_mp$shannon_entropy <- entropy_per_cell

# Summarize entropy per cluster
entropy_summary <- seu_mp@meta.data %>%
    group_by(clusters) %>%
    summarise(
        mean_entropy = mean(shannon_entropy, na.rm = TRUE),
        sd_entropy = sd(shannon_entropy, na.rm = TRUE)
    )
# Plot Shannon entropy
boxplot_entropy <- ggplot(seu_mp@meta.data, aes(x = clusters, y = shannon_entropy, fill = clusters)) +
    geom_boxplot(outlier.size = 0.3) +
    scale_fill_manual(values = sctherapy_colors) +
    theme_bw(base_size = 10) +
    labs(title = "Shannon Entropy per Cluster", y = "Entropy", x = "Cluster") +
    theme_bw(base_size = 9) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = "none"
        )

ggsave("sample_wise/figures/boxplot_entropy.png", plot = boxplot_entropy, width = 5, height = 4, dpi = 300)

pdf("sample_wise/figures/boxplot_entropy.pdf",
    width = 5,
    height = 4
)
boxplot_entropy
dev.off()




# Compute PCA dispersion per cluster
DefaultAssay(seu_mp) <- "RNA"
hvg <- VariableFeatures(seu_mp)
seu_mp <- SketchData(
    object = seu_mp,
    ncells = 50000,
    method = "LeverageScore",
    over.write = TRUE,
    sketched.assay = "sketch_50K",
    features = hvg,
    seed = 120394
)

seu_mp <- ScaleData(seu_mp)
seu_mp <- RunPCA(seu_mp, reduction.name = "pca", npcs = 50)
elbow <- ElbowPlot(seu_mp, ndims = 50)
ggsave("sample_wise/elbowplot.png", plot = elbow)

# Project PCA to the full data
## Set options
options(future.globals.maxSize = 20 * 1024^3)
options(Seurat.object.assay.version = "v5")
seu_mp <- ProjectData(
    object = seu_mp,
    assay = "RNA",
    full.reduction = "pca.full",
    sketched.assay = "sketch_50K",
    sketched.reduction = "pca",
    dims = 1:50
)
# Extract PCA embeddings
DefaultAssay(seu_mp) <- "RNA"
pca_mat <- Embeddings(seu_mp, "pca.full")[, 1:20]

# Add cluster info
meta <- seu_mp@meta.data
meta$cell_id <- rownames(meta)
meta$pca_index <- 1:nrow(meta)

# Combine PCA and cluster info
pca_df <- cbind(meta[, c("clusters", "cell_id")], pca_mat)

# Compute centroid distance per cluster
# Calculate centroids per cluster
centroids <- pca_df %>%
    group_by(clusters) %>%
    summarise(across(starts_with("PC"), mean), .groups = "drop") %>%
    rename_with(~ paste0(., "_centroid"), starts_with("PC"))

# Join centroids back to original cells
pca_cellwise <- pca_df %>%
    left_join(centroids, by = "clusters") %>%
    rowwise() %>%
    mutate(
        distance_to_centroid = sqrt(sum((c_across(starts_with("PC")) - c_across(ends_with("_centroid")))^2))
    ) %>%
    ungroup()

# Plot the PCA dispersion per cluster
boxplot_pca <- ggplot(pca_cellwise, aes(x = clusters, y = distance_to_centroid, fill = clusters)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = sctherapy_colors) +
    labs(
        title = "Distance to Cluster Centroid in PCA Space",
        x = "Cluster",
        y = "Euclidean Distance"
    ) +
    theme_bw(base_size = 9) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = "none"
    )
ggsave("sample_wise/figures/boxplot_pca_dispersion.png", plot = boxplot_pca, width = 5, height = 4, dpi = 300)

pdf("sample_wise/figures/boxplot_pca_dispersion.pdf",
    width = 5,
    height = 4
)
boxplot_pca
dev.off()


# Plot mean transcriptomic and therapeutic heteorgeneity per cluster in a single plot
mean_dispersion <- pca_cellwise %>%
    group_by(clusters) %>%
    summarise(
        mean_distance_to_centroid = mean(distance_to_centroid),
        sd_distance_to_centroid = sd(distance_to_centroid)
    )
colnames(mean_dispersion)[c(2, 3)] <- c("pca_mean", "pca_std")

entropy_df <- seu_mp@meta.data %>%
    group_by(clusters) %>%
    summarise(
        mean_entropy = mean(shannon_entropy, na.rm = TRUE),
        sd_entropy = sd(shannon_entropy, na.rm = TRUE)
    )
colnames(entropy_df)[c(2, 3)] <- c("entropy_mean", "entropy_std")

# Compute a mean Jaccard index per cluster
similarity_matrix <- readRDS("../single_cell/sctherapy/results/jaccard_matrix.rds")
jaccard_dist <- as.dist(1 - similarity_matrix)
cluster_assignment <- readRDS("../single_cell/sctherapy/results/speclustering_reordered.rds")
heterogeneity_by_cluster <- lapply(unique(cluster_assignment), function(cl) {
    samples_in_cl <- names(cluster_assignment[cluster_assignment == cl])

    # Subset distance matrix
    cl_dists <- as.matrix(jaccard_dist)[samples_in_cl, samples_in_cl]

    # Get upper triangle values without diagonal
    dists <- cl_dists[upper.tri(cl_dists)]

    data.frame(
        cluster = cl,
        mean_dist = mean(dists),
        median_dist = median(dists),
        sd_dist = sd(dists),
        n_samples = length(samples_in_cl)
    )
}) %>% bind_rows()

heterogeneity_by_cluster <- heterogeneity_by_cluster %>%
    select(cluster, mean_dist, sd_dist)
colnames(heterogeneity_by_cluster) <- c("clusters", "therapeutic_mean", "therapeutic_std")

joined_df <- mean_dispersion %>%
    left_join(entropy_df, by = "clusters") %>%
    left_join(heterogeneity_by_cluster, by = "clusters")

# Plot the three variables in a single barplot.
# Get min and max per metric
# Get value ranges for scaling
pca_range <- diff(range(joined_df$pca_mean))
entropy_range <- diff(range(joined_df$entropy_mean))
therapeutic_range <- diff(range(joined_df$therapeutic_mean))

plot_df <- joined_df %>%
    mutate(
        # Scale means
        pca_mean_scaled = (pca_mean - min(pca_mean)) / pca_range,
        entropy_mean_scaled = (entropy_mean - min(entropy_mean)) / entropy_range,
        therapeutic_mean_scaled = (therapeutic_mean - min(therapeutic_mean)) / therapeutic_range,

        # Scale stds proportionally to means (preserve ratio)
        pca_std_scaled = pca_std / pca_range,
        entropy_std_scaled = entropy_std / entropy_range,
        therapeutic_std_scaled = therapeutic_std / therapeutic_range
    )


### Plot genomic heterogeneity
### Compute transcriptomic heterogeneity per TC
library(entropy)
cnv_mat <- read.table("../single_cell/cna_metadata/cnv_segments_clones_lvl2_cytobands.tsv")
cnv_mat <- as.matrix(cnv_mat)
cnv_mat <- apply(cnv_mat, c(1,2), as.numeric)
cnv_mat[is.na(cnv_mat)] <- 2

# Discretize (bin) expression to stabilize entropy calculation
# You can use 10 bins (can be tuned)
entropy_per_subclone <- apply(cnv_mat, 2, function(x) {
    entropy(discretize(x, numBins = 5))
})

# Rename subclones
names(entropy_per_subclone) <- gsub("subclone", "", names(entropy_per_subclone))
names(entropy_per_subclone) <- gsub("[^a-zA-Z0-9]", "", names(entropy_per_subclone))

entropy_per_subclone <- as.data.frame(entropy_per_subclone) %>%
    rownames_to_column(var = "subclone_name")
colnames(entropy_per_subclone)[2] <- "entropy"

entropy_df <- seu_mp@meta.data %>%
    select(subclone_name, clusters) %>%
    distinct() %>%
    mutate(subclone_name = gsub("[^a-zA-Z0-9]", "", subclone_name)) %>%
    left_join(entropy_per_subclone, by = "subclone_name")

# Plot Shannon entropy
boxplot_entropy <- ggplot(entropy_df, aes(x = clusters, y = entropy, fill = clusters)) +
    geom_boxplot(outlier.size = 0.3) +
    scale_fill_manual(values = sctherapy_colors) +
    theme_bw(base_size = 10) +
    labs(title = "Genomic heteorogeneity per Cluster (Shanon entropy)", y = "Entropy", x = "Cluster") +
    theme_bw(base_size = 9) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = "none"
        )

ggsave("sample_wise/figures/boxplot_entropy_genomic.png", plot = boxplot_entropy, width = 5, height = 4, dpi = 300)

pdf("sample_wise/figures/boxplot_entropy_genomic.pdf",
    width = 5,
    height = 4
)
boxplot_entropy
dev.off()

# # Plot VlnPlots per cluster
# long_metadata <- long_metadata %>%
#   mutate(
#     Metaprogram = factor(Metaprogram, levels = paste0("MP", 1: 43))
#     )

# vln_boxplot.ucell <- function(data, signature, score){
#     p <- ggplot(data, aes(x = {{signature}}, y = {{score}})) + 
#          geom_violin(aes(fill = {{signature}}), scale = "width") +
#          geom_boxplot(width = 0.1) + theme_bw() +
#          NoLegend()
#     return(p)
# }

# for (i in unique(long_metadata$clusters)){
#     data <- long_metadata[long_metadata$clusters == i,]
#     plot <- vln_boxplot.ucell(data, Metaprogram, UCell_score)
#     ggsave(paste0("sample_wise/vln_boxplot_", i, ".png"), 
#            plot = plot, 
#            width = 15, 
#            height = 8, 
#            dpi = 500)
# }

# vln_boxplot.ucell <- ggplot(long_metadata, aes(x = clusters, y = UCell_score)) + 
#   geom_violin(aes(fill = clusters), scale = "width", trim = FALSE) +
#   geom_boxplot(width = 0.1, outlier.size = 0.5) +
#   facet_wrap(~ Metaprogram, nrow = 5) +  # 5 columns, 2 rows (for 10 clusters)
#   theme_bw(base_size = 9) +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     legend.position = "none"
#   ) +
#   labs(title = "UCell Scores per Metaprogram per Cluster")

# ggsave("sample_wise/combined_vln_boxplot_inverted.png", vln_boxplot.ucell, width = 16, height = 8, dpi = 300)



# This tests if the mean UCell scores for a pathway differ significantly across clusters.

# library(dplyr)
# library(broom)

# anova_results <- long_metadata %>%
#   group_by(Metaprogram) %>%
#   do(tidy(aov(UCell_score ~ clusters, data = .))) %>%
#   filter(term == "clusters") %>%
#   arrange(p.value)

# # View top pathways with most significant variation
# head(anova_results)

# library(purrr)

# # Run Tukey HSD for all MPs and extract significant comparisons
# tukey_results <- long_metadata %>%
#   split(.$Metaprogram) %>%
#   map(~ TukeyHSD(aov(UCell_score ~ clusters, data = .x))$clusters %>%
#         as.data.frame() %>%
#         mutate(comparison = rownames(.), Metaprogram = unique(.x$Metaprogram))) %>%
#   bind_rows() %>%
#   filter(`p adj` < 0.05)  # Only significant differences

# # View top results
# head(tukey_results)


