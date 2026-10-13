library(BPCells)
library(Seurat)
library(dplyr)
library(tidyverse)
library(presto)
library(scCustomize)

setwd("/storage/scratch01/shared/projects/bc-meta/")
source("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")

seu <- readRDS("./single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")

# Leave only malignant cells with clonal information from SCEVAN
full_metadata <- read.table("./functional_nmf/subclone_wise/metadata_subclone_annot.tsv")
metadata_tcca <- read.table("./single_cell/seurat/tcca/tcca_metadata.tsv", sep = "\t", header = TRUE)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/sctherapy/results/")
# Subset the malignant cells with clonal information
seu <- subset(seu, cells = rownames(full_metadata))
seu <- AddMetaData(seu, metadata = full_metadata)

# Add the scTherapy clusters to the metadata
clusters <- readRDS("speclustering_reordered.rds")
clusters_df <- as.data.frame(clusters) %>%
    rownames_to_column("subclone_name")

# Cells without cluster assignation(coming from adrenalnb_rui_chong and
# cell_lines_gabriella_kinker due to the absence of healthy cells for scTherapy
# predictions) will be removed
seu <- subset(seu, subset = subclone_name %in% clusters_df$subclone_name)
seu@meta.data <- seu@meta.data %>%
    rownames_to_column("cell_id") %>%
    left_join(clusters_df, by = "subclone_name") %>%
    column_to_rownames("cell_id")


# Normalize data
seu <- NormalizeData(seu,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)

seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
hvg <- VariableFeatures(seu_before)

# Find differentially expressed features across scTherapy clusters
seu <- JoinLayers(seu)
DefaultAssay(seu) <- "RNA"
Idents(seu) <- "clusters"
seu.markers <- FindAllMarkers(seu, only.pos = TRUE)
write.table(
    as.data.frame(seu.markers),
    "marker_genes/table_markers_stats.tsv",
    sep = "\t", row.names = FALSE
)
saveRDS(seu.markers, "marker_genes/table_markers.rds")

# Compute sketch of 5k cells
seu <- SketchData(
    object = seu,
    ncells = 5000,
    method = "LeverageScore",
    sketched.assay = "sketch",
    features = features_to_include
)

saveRDS(seu, "marker_genes/seu_cluster_markers.rds")


# Save all marker genes per cluster
all_markers <- seu.markers %>%
    select(cluster, gene) %>%
    group_by(cluster) %>%
    mutate(row = row_number()) %>%
    ungroup() %>%
    pivot_wider(names_from = cluster, values_from = gene) %>%
    select(-row) %>%
    as.data.frame()

all_markers <- all_markers[as.character(1:10)]
write.table(top_genes, "marker_genes/all_markers.tsv",
    sep = "\t", row.names = FALSE
)


# Select top 15 genes and create a dotplot
seu.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC >= log2(1.5) & pct.1 >= 0.5 & (pct.1 - pct.2) >= 0.20) %>%
    arrange(desc(avg_log2FC), desc(pct.1)) %>%
    slice_head(n = 15) %>%
    ungroup() -> top15

png("marker_genes/figures/dotplot_top15.png", width = 8, height = 14, res = 300, unit = "in")
p <- scCustomize::Clustered_DotPlot(
    seu,
    features = top15$gene,
    group.by = "clusters",
    colors_use_idents = sctherapy_colors,
    legend_label_size = 11,
    legend_title_size = 11,
    column_label_size = 10,
    row_label_size = 10
)

p + theme(
    plot.margin = unit(c(1, 5, 1, 1), "cm") # top, right, bottom, left
)

dev.off()


# Save top10 for main figure
seu.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC >= log2(1.5) & pct.1 >= 0.5 & (pct.1 - pct.2) >= 0.20) %>%
    arrange(desc(avg_log2FC), desc(pct.1)) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10


pdf("marker_genes/figures/dotplot_top10.pdf", width = 8, height = 12)
p <- scCustomize::Clustered_DotPlot(
    seu,
    features = top10$gene,
    group.by = "clusters",
    colors_use_idents = sctherapy_colors,
    legend_label_size = 11,
    legend_title_size = 11,
    column_label_size = 10,
    row_label_size = 10
)

p + theme(
    plot.margin = unit(c(1, 5, 1, 1), "cm") # top, right, bottom, left
)

dev.off()


# Create a signature with top genes for each cluster
top_genes <- seu.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC >= log2(1.5) & pct.1 >= 0.5 & (pct.1 - pct.2) >= 0.20) %>%
    arrange(desc(avg_log2FC), desc(pct.1)) %>%
    ungroup() %>%
    arrange(cluster)

top_genes.list <- top_genes %>%
    select(cluster, gene) %>%
    mutate(set_name = paste0("Cluster", stringr::str_pad(cluster, width = 2, pad = "0"))) %>%
    group_by(set_name) %>%
    summarise(genes = list(unique(gene)), .groups = "drop") %>%
    deframe()

names(top_genes.list) <- paste0(names(top_genes.list), "_UP")

# Convert to GMT lines
gmt_lines <- lapply(names(top_genes.list), function(sig_name) {
    sig <- paste(c(sig_name, "marker_genes_sctherapy_clusters", top_genes.list[[sig_name]]),
        collapse = "\t"
    )
})

# Write to a single GMT file
writeLines(unlist(gmt_lines), con = "./marker_genes/marker_sigs_filtered.gmt")






# Check TC5 marker genes
seu <- readRDS("./single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
metadata_tcca <- read.table(
    "./single_cell/seurat/tcca/tcca_metadata_h5ad_sctherapy_cells.tsv",
    sep = "\t",
    header = TRUE
)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/sctherapy/results/")
# Subset the malignant cells with clonal information
seu <- subset(seu, cells = metadata_tcca$cell)
metadata_tcca <- metadata_tcca %>% column_to_rownames(var = "cell")
seu <- AddMetaData(seu, metadata = metadata_tcca)

# Subset cells in TC5
seu_tc5 <- subset(seu, subset = therapeutic_cluster == 5)
sort(table(seu_tc5@meta.data$tumor_type), decreasing = TRUE)
x <- seu@meta.data %>%
    filter(therapeutic_cluster == 5) %>%
    select(scevan_subclone, therapeutic_cluster, tumor_type) %>%
    distinct()
sort(table(x$tumor_type), decreasing = TRUE)

# Load TC5 marker gene signature
markers <- readLines("./marker_genes/marker_sigs_filtered.gmt")
markers <- lapply(strsplit(markers, "\t"), function(x) x[-c(1, 2)])
names(markers) <- vapply(
    strsplit(readLines("./marker_genes/marker_sigs_filtered.gmt"), "\t"),
    `[`, character(1), 1
)

# Subset TC5 brain tumors
seu_tc5_brain <- subset(
    seu_tc5,
    subset = refined_tumor_site == "brain" & (tumor_type %in% c("GBM", "LGG", "SKCM"))
)

seu_tc5_brain <- AddModuleScore_UCell(
    seu_tc5_brain,
    features = list(TC5_signature = markers$Cluster05_UP)
)

seu_tc5_brain$group <- case_when(
    seu_tc5_brain$tumor_type == "GBM" ~ "GBM",
    seu_tc5_brain$tumor_type == "LGG" ~ "LGG",
    seu_tc5_brain$tumor_type == "SKCM" ~ "SKCM BrMets",
)

Idents(seu_tc5_brain) <- "group"


vlnplot <- VlnPlot(seu_tc5_brain,
    features = "TC5_signature_UCell",
    pt.size = 0,
    cols = c(
        "GBM" = "#b2509e",
        "LGG" = "#d49dc7",
        "SKCM BrMets" = "#bbd642"
    )
) +
    labs(
        title = "TC5 signature in brain tumors",
        y = "UCell score", x = ""
    ) +
    theme(legend.position = "none")


ggsave("./marker_genes/vln_tc5.pdf", plot = vlnplot, width = 5, height = 6)


# If TC1-3 BRCA and LUAD express more of these genes than BRCA and LUAD in other
# TCs --> marker genes are TC-specific, no cancer type-specific
# Define signatures (replace with your actual genes)
tc_sigs <- list(
    TC1_sig = markers$Cluster01_UP,
    TC2_sig = markers$Cluster02_UP,
    TC3_sig = markers$Cluster03_UP
)

# Subset BRCA and LUAD
seu_sub <- subset(seu, subset = !is.na(therapeutic_cluster) & tumor_type %in% c("BRCA", "LUAD"))
subclone_lvl <- seu@meta.data %>%
    filter(tumor_type %in% c("BRCA", "LUAD")) %>%
    select(scevan_subclone, therapeutic_cluster, tumor_type) %>%
    distinct()
sort(table(subclone_lvl$therapeutic_cluster, subclone_lvl$tumor_type), decreasing = TRUE)

# Compute UCell scores for all signatures at once
seu_sub <- AddModuleScore_UCell(seu_sub, features = tc_sigs)

# Set idents to therapeutic cluster
seu_sub$therapeutic_cluster <- factor(
    paste0("TC", seu_sub$therapeutic_cluster),
    levels = paste0("TC", 1:10)
)
Idents(seu_sub) <- "therapeutic_cluster"

# Plot one violin per signature
plots <- lapply(c("TC1_sig_UCell", "TC2_sig_UCell", "TC3_sig_UCell"), function(sig) {
    VlnPlot(seu_sub,
        features = sig,
        pt.size = 0,
        split.by = "tumor_type",
        cols = c("BRCA" = "#ED2891", "LUAD" = "#D3C3E0")
    ) +
        labs(
            title = gsub("_UCell", "", sig),
            x = "Therapeutic cluster",
            y = "UCell score"
        ) +
        theme(
            legend.position = "right",
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
})

# Arrange plots vertically
vlnplots <- wrap_plots(plots, ncol = 3)
ggsave("./marker_genes/vln_tc1-3.pdf", plot = vlnplots, width = 16, height = 6)
