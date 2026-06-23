library(Seurat)
library(BPCells)
library(UCell)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
#library(openxlsx)
library(matrixStats)
library(fmsb)
library(ggpubr)

setwd("/storage/scratch01/shared/projects/bc-meta/")
source("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")
seu <- readRDS("functional_nmf/sample_wise/seu_mps_sketch.rds")

#### ------------- Transcriptomic Heterogeneity ------------- ####
# 1. Extract variable genes (already calculated by Seurat)
seu$study_sample <- paste0(seu$study, "_", seu$sample)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
hvgs <- VariableFeatures(seu)
logcounts <- t(as.matrix(seu[["RNA"]]$data[hvgs, ])) # cells x genes

# 2. Get a transcriptomic profile per subclone (mean expression of each gene)
subclones <- seu@meta.data$subclone_name
cells_by_subclone <- split(seq_len(nrow(logcounts)), subclones)
profiles <- lapply(cells_by_subclone, function(idx) {
  Matrix::colMeans(logcounts[idx, , drop = FALSE])
})
subclone_profiles <- do.call(cbind, profiles)

# 3. Compute transcriptomic heterogeneity (tITH) per sample taking into account
# or not the therapeutic cluster assignment.
metadata <- read.table(
    "single_cell/seurat/tcca/tcca_metadata_h5ad.tsv", 
    sep = "\t", 
    header = TRUE
    )

subclone_metadata <- metadata %>%
  filter(therapeutic_cluster != "NA") %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  group_by(study, study_sample, scevan_subclone, therapeutic_cluster) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(study_sample) %>%
  mutate(subclone_freq = n_cells / sum(n_cells)) %>%
  ungroup()


compute_tITH <- function(
    subclone_metadata,
    subclone_profiles,
    subclone_col = "scevan_subclone",
    cluster_col = "therapeutic_cluster",
    sample_col = "study_sample",
    freq_col = "subclone_freq",
    unit = c("sample", "TC")
    ) {

    unit  <- match.arg(unit)
    global_dist <- 1 - cor(subclone_profiles, method = "pearson")

    group_col <- if (unit == "sample") sample_col else cluster_col
    groups <- sort(unique(subclone_metadata[[group_col]]))
    groups<- groups[!is.na(groups) & groups != "NA"]

    results <- vapply(groups, function(group) {

        meta_group <- subclone_metadata[subclone_metadata[[group_col]] == group, ]
        subclones  <- meta_group[[subclone_col]]
        subclones  <- subclones[subclones %in% colnames(subclone_profiles)]

        if (length(subclones) < 2) return(0)

        freq <- setNames(meta_group[[freq_col]], meta_group[[subclone_col]])
        freq <- freq[subclones]
        freq <- freq / sum(freq)

        D <- global_dist[subclones, subclones]
        return(as.numeric(freq %*% D %*% freq))

        }, FUN.VALUE = numeric(1))

    names(results) <- groups
    return(results)
}


# Per sample (Fig 3c)
tITH_per_sample <- compute_tITH(subclone_metadata, subclone_profiles, unit = "sample")

# Per TC (Fig 2e)
tITH_per_TC <- compute_tITH(subclone_metadata, subclone_profiles, unit = "TC")


#### ---------------- Genomic Heterogeneity ---------------- #####
# Read CNV matrix per subclone
cnv_mat <- read.table("single_cell/cna_metadata/cnv_segments_clones_lvl2_cytobands.tsv")
colnames(cnv_mat) <- gsub("subclone", "", colnames(cnv_mat))
colnames(cnv_mat) <- gsub("[^a-zA-Z0-9]", "", names(cnv_mat))
cnv_mat <- as.matrix(cnv_mat)
cnv_mat <- apply(cnv_mat, c(1,2), as.numeric)
cnv_mat[is.na(cnv_mat)] <- 2


# Compute distance between subclone
compute_gITH <- function(
    subclone_metadata,
    cnv_mat,
    subclone_col = "scevan_subclone",
    cluster_col = "therapeutic_cluster",
    sample_col = "study_sample",
    freq_col = "subclone_freq",
    unit = c("sample", "TC")
    ) {
  
    cnv_mat <- t(cnv_mat)
    unit <- match.arg(unit)
    L  <- ncol(cnv_mat)
    global_dist <- as.matrix(stats::dist(cnv_mat, method = "euclidean")) / sqrt(L)
    
    # Rename subclones to match cnv mat
    subclone_metadata$scevan_subclone <- gsub("[^a-zA-Z0-9]", "", subclone_metadata$scevan_subclone)
    group_col <- if (unit == "sample") sample_col else cluster_col
    groups <- sort(unique(subclone_metadata[[group_col]]))
    groups <- groups[!is.na(groups) & groups != "NA"]

    results <- vapply(groups, function(group) {

        meta_group <- subclone_metadata[subclone_metadata[[group_col]] == group, ]
        subclones <- meta_group[[subclone_col]]
        subclones <- subclones[subclones %in% rownames(cnv_mat)]

        if (length(subclones) < 2) return(0)

        # Frecuencias directamente desde subclone_metadata
        freq <- setNames(meta_group[[freq_col]], meta_group[[subclone_col]])
        freq <- freq[subclones]
        freq <- freq / sum(freq)  # renormalizar por si algún subclón fue filtrado

        D <- global_dist[subclones, subclones]
        return(as.numeric(freq %*% D %*% freq))

    }, FUN.VALUE = numeric(1))

    names(results) <- groups
    return(results)
}

# Per sample (Fig 3c)
gITH_per_sample <- compute_gITH(subclone_metadata, cnv_mat, unit = "sample")

# Per TC (Fig 2e)
gITH_per_TC <- compute_gITH(subclone_metadata, cnv_mat, unit = "TC")


#### --------------- Therapeutic Heterogeneity per Subclone -------------- ####
# Compute a mean Jaccard index per cluster
similarity_matrix <- readRDS("single_cell/sctherapy/results/jaccard_matrix.rds")
jaccard_dist <- as.matrix(as.dist(1 - similarity_matrix))

compute_thITH <- function(
    subclone_metadata,
    jaccard_dist,
    subclone_col = "scevan_subclone",
    cluster_col = "therapeutic_cluster",
    sample_col = "study_sample",
    freq_col = "subclone_freq",
    unit = c("sample", "TC")
    ) {

  unit <- match.arg(unit)
  group_col <- if (unit == "sample") sample_col else cluster_col
  groups <- sort(unique(subclone_metadata[[group_col]]))
  groups <- groups[!is.na(groups) & groups != "NA"]

  results <- vapply(groups, function(group) {

    meta_group <- subclone_metadata[subclone_metadata[[group_col]] == group, ]
    subclones <- meta_group[[subclone_col]]
    subclones <- subclones[subclones %in% rownames(jaccard_dist)]

    if (length(subclones) < 2) return(0)

    freq <- setNames(meta_group[[freq_col]], meta_group[[subclone_col]])
    freq <- freq[subclones]
    freq <- freq / sum(freq)

    D <- jaccard_dist[subclones, subclones]
    return(as.numeric(freq %*% D %*% freq))

  }, FUN.VALUE = numeric(1))

  names(results) <- groups
  return(results)
}

thITH_per_sample <- compute_thITH(subclone_metadata, jaccard_dist, unit = "sample")
thITH_per_TC <- compute_thITH(subclone_metadata, jaccard_dist, unit = "TC")


# Join the 3 heterogeneities for samples and TCs
sample_lvl_ITH <- data.frame(
    transcriptomic_ITH = tITH_per_sample,
    genomic_ITH = gITH_per_sample,
    therapeutic_ITH = thITH_per_sample
    )
tc_lvl_ITH <- data.frame(
    transcriptomic_ITH = tITH_per_TC,
    genomic_ITH = gITH_per_TC,
    therapeutic_ITH = thITH_per_TC
)

#### ---------------- Boxplots & radar charts per cancer type ---------------- ####
# We add the refined tumor type annotation to ITH
sample_tumor <- metadata %>%
    mutate(study_sample = paste0(study, "_", sample)) %>%
    select(study_sample, tumor_type) %>%
    distinct()

sample_lvl_ITH <-  sample_lvl_ITH %>%
    rownames_to_column(var = "study_sample") %>%
    left_join(sample_tumor, by = "study_sample")

tumor_order <- sample_lvl_ITH %>%
  group_by(tumor_type) %>%
  summarise(median_thITH = median(therapeutic_ITH, na.rm = TRUE)) %>%
  arrange(median_thITH) %>%
  pull(tumor_type)

boxplot_df <- sample_lvl_ITH %>%
  pivot_longer(
    cols      = c(transcriptomic_ITH, genomic_ITH, therapeutic_ITH),
    names_to  = "heterogeneity_type",
    values_to = "heterogeneity_value"
  ) %>%
  mutate(
    tumor_type = factor(tumor_type, levels = tumor_order),
    heterogeneity_type = factor(heterogeneity_type,
      levels = c("therapeutic_ITH", "transcriptomic_ITH", "genomic_ITH"),
      labels = c("Therapeutic", "Transcriptomic", "Genomic")
    )
  )

boxplot1 <- ggplot(
  boxplot_df,
  aes(x = tumor_type, y = heterogeneity_value, 
      fill = tumor_type, color = tumor_type)
) +
  geom_violin(
    alpha = 0.3, 
    color = NA,
    trim  = FALSE
  ) +
  geom_boxplot(
    width         = 0.65,
    outlier.shape = NA,
    linewidth     = 0.3,
    alpha         = 0.6,
    fatten = 1,
    color         = "black"
  ) +
  geom_jitter(
    width  = 0.15,
    size   = 0.8,
    alpha  = 0.6,
    stroke = 0
  ) +
  scale_fill_manual(values  = tumor_type_colors) +
  scale_color_manual(values = tumor_type_colors) +
  facet_wrap(
    ~ heterogeneity_type,
    nrow   = 3,
    scales = "free_y"
  ) +
  labs(x = "Cancer type", y = "Heterogeneity score") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    legend.position  = "none"
  )

ggsave("single_cell/heterogeneity/figures/ith_per_tumor_type.pdf",
  plot   = boxplot1,
  dpi    = 300,
  width  = 6,
  height = 6
)



# To use the fmsb package, we have to add 2 lines to the dataframe: the max and min of
# each variable to show on the plot.
# Format the data for the radarchart.
radar_df <- as.data.frame(t(round(tc_lvl_ITH, digits = 3)))

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

#### ---------------- Boxplots per TC ---------------- ####
# Add TC to sample_lvl_ITH via subclone_metadata
sample_tc <- subclone_metadata %>%
  select(study_sample, therapeutic_cluster) %>%
  distinct()

# Each sample can appear in multiple TCs
sample_lvl_ITH_tc <- sample_lvl_ITH %>%
  left_join(sample_tc, by = "study_sample")

# Orden de TCs por mediana de therapeutic_ITH
tc_order <- sample_lvl_ITH_tc %>%
  group_by(therapeutic_cluster) %>%
  summarise(median_thITH = median(therapeutic_ITH, na.rm = TRUE)) %>%
  arrange(median_thITH) %>%
  pull(therapeutic_cluster)

boxplot_df_tc <- sample_lvl_ITH_tc %>%
  pivot_longer(
    cols = c(therapeutic_ITH, transcriptomic_ITH, genomic_ITH),
    names_to  = "heterogeneity_type",
    values_to = "heterogeneity_value"
  ) %>%
  mutate(
    therapeutic_cluster = factor(therapeutic_cluster, levels = tc_order),
    heterogeneity_type  = factor(heterogeneity_type,
      levels = c("therapeutic_ITH", "transcriptomic_ITH", "genomic_ITH"),
      labels = c("Therapeutic", "Transcriptomic", "Genomic")
    )
  )

boxplot2 <- ggplot(
  boxplot_df_tc,
  aes(x = therapeutic_cluster, y = heterogeneity_value,
      fill = therapeutic_cluster, color = therapeutic_cluster)
) +
  geom_violin(alpha = 0.3, color = NA, trim = FALSE) +
  geom_boxplot(
    width = 0.65, outlier.shape = NA,
    linewidth = 0.3, alpha = 0.6,
    fatten = 1, color = "black"
  ) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.6, stroke = 0) +
  scale_fill_manual(values  = sctherapy_colors) +
  scale_color_manual(values = sctherapy_colors) +
  facet_wrap(
    ~ heterogeneity_type,
    ncol = 3,
    scales = "free_y"
  ) +
  labs(x = "Therapeutic cluster", y = "Heterogeneity score") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x      = element_text(angle = 0),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", hjust = 0),
    legend.position  = "none"
  )

ggsave("single_cell/heterogeneity/figures/ith_per_TC.pdf",
  plot   = boxplot2,
  dpi    = 300,
  width  = 8,
  height = 3
)


# Add number of cells per sample
seu <- readRDS("functional_nmf/sample_wise/seu_mps_sketch.rds")
num_cells_sample <- seu@meta.data %>%
    mutate(study_sample = paste0(study, "_", sample)) %>%
    group_by(study_sample) %>%
    summarise(n_cells = n())
corr_genomic_therapeutic <- sample_lvl_ITH  %>%
    inner_join(num_cells_sample, by = "study_sample")

# Add number of drugs per sample
sctherapy <- read.table(
    "single_cell/sctherapy/results/drug_response_subclone_sctherapy.tsv",
    sep = "\t", header = TRUE
)
num_drugs_sample <- sctherapy %>%
    mutate(study_sample = paste0(Study, "_", Sample)) %>%
    group_by(study_sample) %>%
    summarise(n_drugs = n_distinct(Drug.Name))
corr_genomic_therapeutic <- corr_genomic_therapeutic %>%
    inner_join(num_drugs_sample, by = "study_sample")


# Join clinical data to genomic and therapeutic ITH table
clinical_data <- metadata %>%
    mutate(study_sample = paste0(study, "_", sample)) %>%
    select(
        study_sample,
        #tumor_type,
        sample_type,
        treated
    ) %>%
    distinct()
    
corr_genomic_therapeutic <- corr_genomic_therapeutic %>%
    left_join(clinical_data, by = "study_sample")

# Plot distribution of number of drugs per sample
plot <- ggplot(corr_genomic_therapeutic, aes(x = n_drugs)) +
    geom_histogram(binwidth = 5, fill = "#3182bd", color = "black", alpha = 0.7) +
    geom_vline(xintercept = c(20, 40), linetype = "dashed", color = "red") +
    labs(
        title = "Distribución del número de fármacos efectivos por muestra",
        x = "Número de fármacos (n_drugs)",
        y = "Frecuencia"
    ) +
    theme_bw(base_size = 14)
ggsave(
    "single_cell/genomic_ith/figures/num_drugs_distribution.png",
    plot = plot,
    width = 15,
    height = 6,
    dpi = 300
)

# Plot correlation between genomic and therapeutic heterogeneity
theme <- theme_minimal() +
    theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", face = "bold"),
        legend.title = element_text(color = "black", face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )

corr_plot_faceted <- corr_genomic_therapeutic %>%
    filter(sample_type %in% c("p", "m")) %>%
    mutate(
        sample_type = factor(
            sample_type,
            levels = c("p", "m"),
            labels = c("Primary", "Metastasis")
        ),
        n_drugs_cat = cut(
            n_drugs,
            breaks = c(0, 30, 40, 50, 60, Inf),
            labels = c("<30", "30-40", "40-50", "50-60", ">60"),
            right = TRUE
        )
    ) %>%
    ggplot(aes(
        x = therapeutic_ITH,
        y = genomic_ITH,
        color = tumor_type,
        size = n_drugs_cat
    )) +
    geom_point(alpha = 0.8) +
    facet_wrap(~sample_type, nrow = 2) +
    labs(
        x = "Therapeutic ITH",
        y = "Genomic ITH",
        color = "Cancer type",
        size = "Number of effective drugs"
    ) +
    scale_color_manual(values = tumor_type_colors) +
    scale_size_manual(
        values = c(2, 3, 4.5, 6, 8), # tamaños crecientes por categoría
        drop = FALSE
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 9)) +
    theme_bw(base_size = 14) +
    theme_bw(base_size = 14) +
    theme(
        panel.border = element_blank(),
        axis.line = element_line(
            color = "black",
            linewidth = 0.6
        ),
        panel.grid = element_blank(),

        strip.background = element_blank(),
        strip.text = element_text(
            face = "bold",
            size = 14,
            hjust = 0
        ),

        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.background = element_blank(),
        legend.box.background = element_blank()
    )
ggsave(
    "single_cell/genomic_ith/figures/genomic_vs_therapeutic_ith_primary_vs_met.pdf",
    plot = corr_plot_faceted,
    width = 8,
    height = 10,
    dpi = 300
)


# Interactive version of the plot
library(ggplot2)
library(plotly)
library(dplyr)

corr_plot_faceted_interactive <- corr_genomic_therapeutic %>%
    filter(sample_type %in% c("p", "m")) %>%
    mutate(
        sample_type = factor(
            sample_type,
            levels = c("p", "m"),
            labels = c("Primary", "Metastasis")
        ),
        n_drugs_cat = cut(
            n_drugs,
            breaks = c(0, 30, 40, 50, 60, Inf),
            labels = c("<30", "30-40", "40-50", "50-60", ">60"),
            right = TRUE
        )
    ) %>%
    ggplot(aes(
        x = therapeutic_ITH,
        y = genomic_ITH,
        color = tumor_type,
        size = n_drugs_cat,
        text = paste0(
            "Sample: ", study_sample, "<br>",
            "Cancer type: ", tumor_type, "<br>",
            "Genomic ITH: ", round(genomic_ITH, 2), "<br>",
            "Therapeutic ITH: ", round(therapeutic_ITH, 2), "<br>",
            "Effective drugs: ", n_drugs
        )
    )) +
    geom_point(alpha = 0.8) +
    facet_wrap(~sample_type, ncol = 2) +
    labs(
        x = "Therapeutic ITH",
        y = "Genomic ITH"
    ) +
    scale_color_manual(values = tumor_type_colors) +
    scale_size_manual(
        values = c(2, 3, 4.5, 6, 8),
        drop = FALSE
    ) +
    guides(color = "none", size = "none") + # Esto quita las leyendas
    theme_bw(base_size = 14) +
    theme(
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
        strip.text = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid = element_blank(),
        legend.position = "none" # También puedes usar esto
    )

# Convertir a plotly
corr_plot_faceted_interactive <- ggplotly(corr_plot_faceted_interactive, tooltip = "text")

htmlwidgets::saveWidget(
    corr_plot_faceted_interactive,
    "single_cell/genomic_ith/figures/corr_plot_faceted_interactive.html",
    selfcontained = TRUE
)


# Compute correlation matrix
# Transcriptomic vs Therapeutic
cor.test(
    sample_lvl_ITH$transcriptomic_ITH,
    sample_lvl_ITH$therapeutic_ITH,
    method = "pearson"
)

# Transcriptomic vs Genomic
cor.test(
    sample_lvl_ITH$transcriptomic_ITH,
    sample_lvl_ITH$genomic_ITH,
    method = "pearson"
)

# Therapeutic vs Genomic
cor.test(
    sample_lvl_ITH$therapeutic_ITH,
    sample_lvl_ITH$genomic_ITH,
    method = "pearson"
)