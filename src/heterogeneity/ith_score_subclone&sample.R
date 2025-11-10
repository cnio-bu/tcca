library(Seurat)
library(BPCells)
library(UCell)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(openxlsx)
library(entropy)
library(matrixStats)
library(fmsb)
library(ggpubr)

setwd("/storage/scratch01/shared/projects/bc-meta/")
source("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")
seu <- readRDS("functional_nmf/sample_wise/seu_mps_sketch.rds")

#### ------------- Transcriptomic Heterogeneity per Subclone ------------- ####
# Get log-normalized expression matrix
expr_mat <- as.matrix(seu[["RNA"]]$data)

# 1. Extract variable genes (already calculated by Seurat)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
variable_genes <- VariableFeatures(seu)
expr_var <- as.matrix(seu[["RNA"]]$data[variable_genes, ])

# 2. Pre-calculate cell masks per subclone
subclones <- unique(seu@meta.data$subclone_name)
cell_masks <- lapply(subclones, function(sub) {
    seu@meta.data$subclone_name == sub
})
names(cell_masks) <- subclones

# 3. Calculate heterogeneity
heterogeneity_results <- data.frame(
    subclone = subclones,
    mean_cv = NA,
    median_cv = NA,
    n_cells = NA,
    n_variable_genes = NA
)

for (i in seq_along(subclones)) {
    sub <- subclones[i]
    cells <- cell_masks[[sub]]
    n_cells <- sum(cells)

    if (n_cells < 3) { # Minimum 3 cells for reliable CV
        heterogeneity_results$n_cells[i] <- n_cells
        next
    }

    expr_sub <- expr_var[, cells, drop = FALSE]

    # Calculate CV per gene
    gene_means <- rowMeans(expr_sub)
    gene_sds <- rowSds(expr_sub)

    # CV = sd/mean, only for expressed genes
    expressed <- gene_means > 0.01 # Minimum threshold
    cv_genes <- gene_sds[expressed] / gene_means[expressed]

    # Filter extreme outliers (optional but recommended)
    cv_genes <- cv_genes[cv_genes < quantile(cv_genes, 0.99, na.rm = TRUE)]

    heterogeneity_results$mean_cv[i] <- mean(cv_genes, na.rm = TRUE)
    heterogeneity_results$median_cv[i] <- median(cv_genes, na.rm = TRUE)
    heterogeneity_results$n_cells[i] <- n_cells
    heterogeneity_results$n_variable_genes[i] <- sum(expressed)
}

# write.table(heterogeneity_results,
#     file = "sample_wise/single_cell/heterogeneity/transcriptomic_heterogeneity_per_subclone.tsv",
#     sep = "\t",
#     quote = FALSE,
#     row.names = FALSE
# )
# 4. Add metadata and visualize
seu@meta.data$transcriptomic_ith_cv <- heterogeneity_results$mean_cv[
    match(seu@meta.data$subclone_name, heterogeneity_results$subclone)
]


#### --------------- Therapeutic Heterogeneity per Subclone -------------- ####
# Compute a mean Jaccard index per cluster
similarity_matrix <- readRDS("single_cell/sctherapy/results/jaccard_matrix.rds")
jaccard_dist <- as.dist(1 - similarity_matrix)
cluster_assignment <- readRDS("single_cell/sctherapy/results/speclustering_reordered.rds")

cluster_list <- split(names(cluster_assignment), cluster_assignment)
subclone_heterogeneity <- sapply(names(cluster_assignment), function(s) {
    cl <- cluster_assignment[s]
    cluster_samples <- names(cluster_assignment[cluster_assignment == cl])
    cluster_samples <- setdiff(cluster_samples, s) # exclude self
    sample_distance <- mean(1 - similarity_matrix[s, cluster_samples]) # mean distance
    return(sample_distance)
})

# Add to seurat object
therapeutic_het <- enframe(subclone_heterogeneity, name = "Subclone", value = "therapeutic_heterogeneity")
seu@meta.data <- seu@meta.data %>%
    left_join(therapeutic_het, by = c("subclone_name" = "Subclone"))





#### ---------------- Genomic Heterogeneity per Subclone ---------------- ####
# Read CNV matrix per subclone
library(entropy)
cnv_mat <- read.table("single_cell/cna_metadata/cnv_segments_clones_lvl2_cytobands.tsv")
cnv_mat <- as.matrix(cnv_mat)
cnv_mat <- apply(cnv_mat, c(1,2), as.numeric)
cnv_mat[is.na(cnv_mat)] <- 2


# Discretize (bin) expression to stabilize entropy calculation
# You can use 10 bins (can be tuned)
genomic_entropy <- apply(cnv_mat, 2, function(x) {
    entropy(discretize(x, numBins = 5))
})

# Rename subclones
names(genomic_entropy) <- gsub("subclone", "", names(genomic_entropy))
names(genomic_entropy) <- gsub("[^a-zA-Z0-9]", "", names(genomic_entropy))

genomic_entropy <- enframe(genomic_entropy,
    name = "Subclone", 
    value = "shannon_entropy_genomic"
)

seu@meta.data <- seu@meta.data %>%
    mutate(subclone_match = gsub("[^a-zA-Z0-9]", "", subclone_name)) %>%
    left_join(genomic_entropy, by = c("subclone_match" = "Subclone")) %>%
    select(-subclone_match)


# Save subclone-level heterogeneity scores
ith_subclone <- seu@meta.data %>%
    select(
        subclone_name,
        transcriptomic_ith_cv,
        shannon_entropy_genomic,
        therapeutic_heterogeneity
    ) %>%
    distinct()
write.table(ith_subclone,
    file = "single_cell/heterogeneity/ith_scores_per_subclone.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)



#### ---------------- Boxplots & radar charts per TC ---------------- ####
# We add the refined tumor type annotation to the metadata
metadata <- read.table("single_cell/seurat/v5/tcca_metadata.tsv",
    sep = "\t",
    header = TRUE
)
seu@meta.data <- seu@meta.data %>%
    left_join(select(metadata, cell, refined_tumor_type), by = "cell")

subclone_heterogeneity <- seu@meta.data %>%
    group_by(subclone_name, clusters, refined_tumor_type) %>%
    summarize(
        transcriptomic_mean = mean(transcriptomic_ith_cv, na.rm = TRUE),
        genomic_mean = mean(shannon_entropy_genomic, na.rm = TRUE),
        therapeutic_mean = mean(therapeutic_heterogeneity, na.rm = TRUE)
    )

boxplot_df <- subclone_heterogeneity %>%
    pivot_longer(
        cols = c(transcriptomic_mean, genomic_mean, therapeutic_mean),
        names_to = "heterogeneity_type", values_to = "heterogeneity_value"
    )

th <- theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "black")
)

boxplot_df$heterogeneity_type <- factor(
    boxplot_df$heterogeneity_type,
    levels = c("therapeutic_mean", "transcriptomic_mean", "genomic_mean")
)
tc_order <- seu@meta.data %>%
    select(clusters, subclone_name, therapeutic_heterogeneity) %>%
    distinct() %>%
    mutate(clusters = fct_reorder(
        clusters,
        therapeutic_heterogeneity,
        .fun = median,
        .desc = FALSE
    )) %>%
    pull(clusters)
boxplot_df$clusters <- factor(boxplot_df$clusters, levels = levels(tc_order))
boxplot1 <- ggplot(boxplot_df, aes(x = clusters, y = heterogeneity_value, fill = clusters)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.7) +
    labs(
        x = "scTherapy cluster",
        y = "Heterogeneity score",
        fill = "Heterogeneity"
    ) +
    facet_wrap(~heterogeneity_type, ncol = 3, scales = "free_y", labeller = as_labeller(c(
                transcriptomic_mean = "Transcriptomic",
                genomic_mean = "Genomic",
                therapeutic_mean = "Therapeutic"))) +
    scale_fill_manual(values = sctherapy_colors) +
    theme_bw() +
    th +
    theme(legend.position = "none")

ggsave("single_cell/heterogeneity/figures/ith_per_tc.pdf",
    plot = boxplot1,
    dpi = 300,
    width = 12, height = 4
)




# To use the fmsb package, we have to add 2 lines to the dataframe: the max and min of
# each variable to show on the plot.
# Format the data for the radarchart.
radar_df <- subclone_heterogeneity %>%
    select(
        clusters, transcriptomic_mean,
        therapeutic_mean,
        genomic_mean
    ) %>%
    group_by(clusters) %>%
    summarize(
        transcriptomic = mean(transcriptomic_mean, na.rm = TRUE),
        genomic = mean(genomic_mean, na.rm = TRUE),
        therapeutic = mean(therapeutic_mean, na.rm = TRUE)
    ) %>%
    # mutate(
    #     transcriptomic = scale(transcriptomic)[, 1],
    #     genomic = scale(genomic)[, 1],
    #     therapeutic = scale(therapeutic)[, 1]
    # ) %>%
    column_to_rownames(var = "clusters") %>%
    t() %>%
    as.data.frame()


radar_df <- round(radar_df, digits = 3)

# Plot each row (heterogeneity type) as its own radar chart
custom_titles <- c("Transcriptomic", "Genomic", "Therapeutic")
colors_border <- c("#d96157", "#447dac", "#53744c")
colors_in <- alpha(colors_border, 0.3)

# Set up PDF or PNG
pdf("single_cell/heterogeneity/figures/radarchart_ith_per_tc.pdf",
    width = 10, height = 4
)
par(mfrow = c(1, 3))
par(mar = c(2, 1, 3, 1))
par(oma = c(0, 0, 2, 0))
for (i in seq_along(custom_titles)) {
    df_i <- radar_df[i, , drop = FALSE]
    df_i <- rbind(rep(max(df_i), ncol(df_i)), rep(min(df_i), ncol(df_i)), df_i)

    radarchart(df_i,
        vlabels = colnames(df_i),
        pcol = colors_border[i],
        pfcol = colors_in[i],
        plwd = 4,
        plty = 1,
        cglcol = "grey", cglty = 1, axislabcol = "black", cglwd = 0.8,
        vlcex = 1.2,
        title = custom_titles[i]
    )
}
dev.off()


#### ------------ Boxplots & radar charts per tumor type -------------- ####
cancer_order <- seu@meta.data %>%
    select(clusters, refined_tumor_type, subclone_name, therapeutic_heterogeneity) %>%
    distinct() %>%
    mutate(refined_tumor_type = fct_reorder(
        refined_tumor_type,
        therapeutic_heterogeneity,
        .fun = median,
        .desc = FALSE
    )) %>%
    pull(refined_tumor_type)

boxplot_df$refined_tumor_type <- factor(boxplot_df$refined_tumor_type, levels = levels(cancer_order))

boxplot2 <- ggplot(boxplot_df, aes(
    x = factor(refined_tumor_type),
    y = heterogeneity_value,
    fill = factor(refined_tumor_type)
)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.7) +
    labs(
        x = "Tumor type",
        y = "Heterogeneity score",
        fill = "Heterogeneity"
    ) +
    facet_wrap(~heterogeneity_type, ncol = 1, scales = "free_y", labeller = as_labeller(c(
        transcriptomic_mean = "Transcriptomic",
        genomic_mean = "Genomic",
        therapeutic_mean = "Therapeutic"
    ))) +
    scale_fill_manual(values = tumor_type_colors) +
    theme_bw() +
    th +
    theme(legend.position = "none")



ggsave("single_cell/heterogeneity/figures/ith_per_tumor_type_facet.pdf",
    plot = boxplot2,
    dpi = 300,
    width = 12, height = 10
)

# Radar chart per cancer type
# Format the data for the radarchart.
radar_df <- subclone_heterogeneity %>%
    select(
        refined_tumor_type, transcriptomic_mean,
        therapeutic_mean,
        genomic_mean
    ) %>%
    group_by(refined_tumor_type) %>%
    summarize(
        transcriptomic = median(transcriptomic_mean, na.rm = TRUE),
        genomic = median(genomic_mean, na.rm = TRUE),
        therapeutic = median(therapeutic_mean, na.rm = TRUE)
    ) %>%
    # mutate(
    #     transcriptomic = scale(transcriptomic)[, 1],
    #     genomic = scale(genomic)[, 1],
    #     therapeutic = scale(therapeutic)[, 1]
    # ) %>%
    column_to_rownames(var = "refined_tumor_type") %>%
    t() %>%
    as.data.frame()


radar_df <- round(radar_df, digits = 3)

# Plot each row (heterogeneity type) as its own radar chart
custom_titles <- c("Transcriptomic", "Genomic", "Therapeutic")
colors_border <- c("#d96157", "#447dac", "#53744c")
colors_in <- alpha(colors_border, 0.3)

# Set up PDF or PNG
pdf("single_cell/heterogeneity/figures/radarchart_ith_per_cancertype.pdf",
    width = 10, height = 4
)
par(mfrow = c(1, 3))
par(mar = c(2, 1, 3, 1))
par(oma = c(0, 0, 2, 0))
for (i in seq_along(custom_titles)) {
    df_i <- radar_df[i, , drop = FALSE]
    df_i <- rbind(rep(max(df_i), ncol(df_i)), rep(min(df_i), ncol(df_i)), df_i)

    radarchart(df_i,
        vlabels = colnames(df_i),
        pcol = colors_border[i],
        pfcol = colors_in[i],
        plwd = 4,
        plty = 1,
        cglcol = "grey", cglty = 1, axislabcol = "black", cglwd = 0.8,
        vlcex = 1.2, # 🔹 bigger cluster labels
        title = custom_titles[i] # 🔹 custom plot title
    )
}
dev.off()




#### ---------------- Correlations between heterogeneities ---------------- ####
heterogeneity_df <- seu@meta.data %>%
    select(
        subclone_name,
        transcriptomic_ith_cv,
        therapeutic_heterogeneity,
        shannon_entropy_genomic
    ) %>%
    distinct()

# Generate all pairwise combinations (without repetition)
var_pairs <- combn(names(heterogeneity_df)[2:4], 2, simplify = FALSE)

# Create a list of ggscatter plots
plots <- map(var_pairs, function(vars) {
    xvar <- vars[1]
    yvar <- vars[2]
    ggscatter(
        heterogeneity_df,
        x = xvar,
        y = yvar,
        add = "reg.line",
        conf.int = TRUE,
        cor.coef = TRUE,
        cor.method = "spearman",
        cor.coef.size = 5
    ) +
        ggtitle(paste(xvar, "vs", yvar))
})

# Combine all scatter plots into a single figure
combined_plot <- ggpubr::ggarrange(
    plotlist = plots,
    ncol = 3, nrow = 1,
    labels = LETTERS[1:length(plots)]
)

# Display
pdf("single_cell/heterogeneity/figures/heterogeneity_correlations.pdf",
    width = 14,
    height = 4
)
combined_plot
dev.off()


# Compute correlation matrix
# Transcriptomic vs Therapeutic
cor.test(
    heterogeneity_df$transcriptomic_ith_cv,
    heterogeneity_df$therapeutic_heterogeneity,
    method = "pearson"
)

# Transcriptomic vs Genomic
cor.test(
    heterogeneity_df$transcriptomic_ith_cv,
    heterogeneity_df$shannon_entropy_genomic,
    method = "pearson"
)

# Therapeutic vs Genomic
cor.test(
    heterogeneity_df$therapeutic_heterogeneity,
    heterogeneity_df$shannon_entropy_genomic,
    method = "pearson"
)





#### -------------------------- ITH per sample --------------------------- ####
# Compute frequency of subclones per sample
colnames(similarity_matrix) <- gsub("[^a-zA-Z0-9]", "", colnames(similarity_matrix))
rownames(similarity_matrix) <- gsub("[^a-zA-Z0-9]", "", rownames(similarity_matrix))

subclone_freq <- metadata %>%
    mutate(scevan_subclone = gsub("[^a-zA-Z0-9]", "", scevan_subclone)) %>%
    filter(malignancy == "True" & scevan_subclone %in% colnames(similarity_matrix)) %>%
    group_by(study_sample, scevan_subclone) %>%
    summarise(n_cells = n()) %>%
    group_by(study_sample) %>%
    mutate(p = n_cells / sum(n_cells)) %>%
    ungroup()

# Compute transcriptomic heterogeneity per sample (weighted by number of cells in each subclone)
transcriptomic_ith <- heterogeneity_results %>%
    left_join(
        seu@meta.data %>%
            mutate(study_sample = paste0(study,"_", sample)) %>%
            select(subclone_name, study_sample) %>%
            distinct(),
        by = c("subclone" = "subclone_name")
    ) %>%
    filter(!is.na(mean_cv), n_cells >= 3) %>%
    group_by(study_sample) %>%
    summarise(
        transcriptomic_ith = weighted.mean(mean_cv, w = n_cells, na.rm = TRUE),
        n_subclones = n(),
        total_cells = sum(n_cells)
    )


# Compute genomic heterogeneity per sample
colnames(cnv_mat) <- gsub("subclone", "", gsub("[^a-zA-Z0-9]", "", colnames(cnv_mat)))

samples <- unique(subclone_freq$study_sample)
entropy_sample <- setNames(numeric(length(samples)), samples)

prop_cnv_by_subclon <- lapply(colnames(cnv_mat), function(subclone) {
    cnv_subclone <- cnv_mat[, subclone]
    bins <- discretize(cnv_subclone, numBins = 5)
    prop <- bins / sum(bins)
    return(prop)
})

names(prop_cnv_by_subclon) <- colnames(cnv_mat)

for (s in samples) {
    subclones <- subclone_freq$scevan_subclone[subclone_freq$study_sample == s]
    fracs <- subclone_freq$p[subclone_freq$study_sample == s]
    mat <- do.call(rbind, prop_cnv_by_subclon[subclones])
    weighted_prop <- colSums(mat * fracs)
    entropy_sample[s] <- entropy(weighted_prop, unit = "log2")
}

genomic_ith <- data.frame(
    study_sample = names(entropy_sample),
    genomic_ith = entropy_sample,
    row.names = NULL,
    stringsAsFactors = FALSE
)

# Compute therapeutic heterogeneity per sample (Rao's Quadratic Entropy, Rao Q)
dist_mat <- 1 - similarity_matrix

compute_rao <- function(p_vec, dist_mat) {
    subclones <- names(p_vec)
    dist_sub <- dist_mat[subclones, subclones, drop = FALSE]
    p <- matrix(p_vec, ncol = 1)
    as.numeric(t(p) %*% dist_sub %*% p)
}

therapeutic_ith <- subclone_freq %>%
    group_by(study_sample) %>%
    summarise(
        rao_q = {
            p_vec <- setNames(p, scevan_subclone)
            compute_rao(p_vec, dist_mat)
        }
    )


# Join the 3 heterogeneities
heterogeneity_per_sample <- transcriptomic_ith %>%
    select(study_sample, transcriptomic_ith) %>%
    left_join(genomic_ith, by = "study_sample") %>%
    left_join(therapeutic_ith, by = "study_sample")
# Save table with ITH per sample
write.table(heterogeneity_per_sample,
    file = "single_cell/heterogeneity/ith_scores_per_sample.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# Boxplots of sample ITH per cancer type
heterogeneity_per_sample <- heterogeneity_per_sample %>%
    left_join(
        seu@meta.data %>%
            mutate(study_sample = paste0(study, "_", sample)) %>%
            select(study_sample, refined_tumor_type) %>%
            distinct(),
        by = "study_sample"
    )

boxplot_df <- heterogeneity_per_sample %>%
    pivot_longer(
        cols = c(transcriptomic_ith, genomic_ith, rao_q),
        names_to = "heterogeneity_type", values_to = "heterogeneity_value"
    )

th <- theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "black")
)

boxplot_df$heterogeneity_type <- factor(
    boxplot_df$heterogeneity_type,
    levels = c("rao_q", "transcriptomic_ith", "genomic_ith")
)
cancer_order <- heterogeneity_per_sample %>%
    select(refined_tumor_type, study_sample, rao_q) %>%
    distinct() %>%
    mutate(refined_tumor_type = fct_reorder(
        refined_tumor_type,
        rao_q,
        .fun = median,
        .desc = FALSE
    )) %>%
    pull(refined_tumor_type)

boxplot_df$refined_tumor_type <- factor(boxplot_df$refined_tumor_type, levels = levels(cancer_order))

boxplot2 <- ggplot(boxplot_df, aes(
    x = factor(refined_tumor_type),
    y = heterogeneity_value,
    fill = factor(refined_tumor_type)
)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.7) +
    labs(
        x = "Tumor type",
        y = "Heterogeneity score",
        fill = "Heterogeneity"
    ) +
    facet_wrap(~heterogeneity_type, ncol = 1, scales = "free_y", labeller = as_labeller(c(
        transcriptomic_ith = "Transcriptomic",
        genomic_ith = "Genomic",
        rao_q = "Therapeutic"
    ))) +
    scale_fill_manual(values = tumor_type_colors) +
    theme_bw() +
    th +
    theme(legend.position = "none")

ggsave("single_cell/heterogeneity/figures/sample_ith_per_cancertype.pdf",
    plot = boxplot2,
    dpi = 300,
    width = 12, height = 10
)


# Correlation between transcriptomic, genomic and therapeutic heterogeneity per sample
# Generate all pairwise combinations (without repetition)
var_pairs <- combn(colnames(heterogeneity_per_sample)[2:4], 2, simplify = FALSE)

# Create a list of ggscatter plots
plots <- map(var_pairs, function(vars) {
    xvar <- vars[1]
    yvar <- vars[2]
    ggscatter(
        heterogeneity_per_sample,
        x = xvar,
        y = yvar,
        add = "reg.line",
        conf.int = TRUE,
        cor.coef = TRUE,
        cor.method = "spearman",
        cor.coef.size = 5
    ) +
        ggtitle(paste(xvar, "vs", yvar))
})

# Combine all scatter plots into a single figure
combined_plot <- ggpubr::ggarrange(
    plotlist = plots,
    ncol = 3, nrow = 1,
    labels = LETTERS[1:length(plots)]
)

# Display
pdf("single_cell/heterogeneity/figures/ith_correlations_sample_wise.pdf",
    width = 14,
    height = 4
)
combined_plot
dev.off()
