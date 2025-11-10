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

# Load seurat object with malignant cells
seu_lvl2 <- readRDS("../single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
malignant <- subset(seu_lvl2, subset = malignancy == TRUE)

# Load the list of metaprograms
mp_list <- readRDS("sample_wise/metaprograms_cpm/mp_list_reordered.rds")
names(mp_list) <- paste0("MP", 1:length(mp_list))
mp_ucell <- AddModuleScore_UCell(malignant, features = mp_list)

# Save seurat object with UCell scores in metadata
saveRDS(mp_ucell, "sample_wise/seurat_mps_ucell.rds")
# Save table of UCell scores
write.table(mp_ucell@meta.data[, paste0("MP", 1: 43, "_UCell")], 
            "sample_wise/mps_ucell_scores.tsv")

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

# Compute mean MPs enrichment per subclone
seu_mp <- readRDS("sample_wise/seurat_mps_ucell.rds")
tcca_metadata <- read.table(
    "../single_cell/seurat/v5/tcca_metadata.tsv",
    sep = "\t",
    header = TRUE
)
mps <- seu_mp@meta.data %>%
    select(cell, contains("UCell")) %>%
    select(-CIN70_UCell)

mp_cols <- grep("_UCell$", colnames(mps), value = TRUE)

# Assign top MP per subclone based on mean UCell scores
temp <- tcca_metadata %>%
    left_join(mps, by = "cell") %>%
    filter(malignancy == "True")

subclone_mps <- temp %>%
    filter(scevan_subclone != "") %>%
    group_by(scevan_subclone) %>%
    summarise(across(all_of(mp_cols), mean, na.rm = TRUE))

subclone_mps <- subclone_scores %>%
    rowwise() %>%
    mutate(
        top_MP = mp_cols[which.max(c_across(all_of(mp_cols)))]
    ) %>%
    ungroup() %>%
    as.data.frame()

# Translate top_MP to MP1, MP2, etc
mp_names_clean <- c(
    "CellCycle.G2M",
    "CellCycle.G1S",
    "CellCycle.HMG-rich",
    "CellCycle.ChromatinRemodeling",
    "CellCycle..Chromatin",
    "CellCycle.DNARepair",
    "Oncogenic.MYC",
    "Stress.IEGs.AP1",
    "Stress.ISR",
    "Stress.Hypoxia",
    "Stress.Metabolic",
    "Stress.Detoxification",
    "Stress.OxidativeStress",
    "Inflammation.Interferon.MHCII",
    "Inflammation.ReactiveEpithelia",
    "Inflammation.TNFA.NFkB",
    "EMT.partialEMT",
    "EMT.EMT_I",
    "EMT.EMT_II",
    "EMT.Mesenchymal-like",
    "EMT.EMT_III",
    "CellularPlasticity.ActiveSignaling",
    "CellularPlasticity.Post-transcriptionalChromatin",
    "CellularPlasticity.EpithelialRemodeling",
    "ProteinRegulation.ProteasomeDegradation",
    "ProteinRegulation.ProteinMaturation",
    "ProteinRegulation.UPR",
    "ProteinRegulation.ProteinTranslation",
    "EpithelialSenescence",
    "MitochondrialRespiration",
    "Cilia",
    "LineageSpecific.Hemato.InflammatoryMyeloid",
    "LineageSpecific.Hemato.HPSCs",
    "LineageSpecific.Hemato.Neutrophil",
    "LineageSpecific.Hemato.APC-MHCII",
    "LineageSpecific.Hemato.MastCells",
    "LineageSpecific.Hemato.RBCs",
    "LineageSpecific.Neural.Astrocytes",
    "LineageSpecific.Neural.OPCs-NPCs",
    "LineageSpecific.Melanocyte-Pigmentation",
    "LineageSpecific.Prostate-Secretory",
    "LineageSpecific.Urothelial-Secretory",
    "LineageSpecific.Gastrointestinal-Secretory"
)
names(mp_names_clean) <- mp_cols

subclone_mps$top_MP_clean <- mp_names_clean[subclone_mps$top_MP]
subclone_mps <- subclone_mps %>%
    select(scevan_subclone, top_MP, top_MP_clean)
write.table(subclone_mps, "sample_wise/subclone_top_mps.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# Assign cells to MP with highest score
malignant_mp <- tcca_metadata %>%
    left_join(mps, by = "cell") %>%
    filter(malignancy == "True") %>%
    rowwise() %>%
        mutate(
            top_MP = mp_cols[which.max(c_across(all_of(mp_cols)))]
        ) %>%
        ungroup() %>%
        as.data.frame()

malignant_mp$top_MP_clean <- mp_names_clean[malignant_mp$top_MP]

malignant_mp <- malignant_mp %>%
    mutate(top_MP_clean = str_c(str_remove(top_MP, "_UCell$"), "_", top_MP_clean))

tcca_metadata <- tcca_metadata %>%
    left_join(
        malignant_mp %>%
            select(cell, top_MP_clean),
        by = "cell"
    )
