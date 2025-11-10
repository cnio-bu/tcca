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
seu <-  ScaleData(seu, features = rownames(seu))
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
gmt_lines <- lapply(names(top_genes.list), function(sig_name){
    sig <- paste(c(sig_name, "marker_genes_sctherapy_clusters", top_genes.list[[sig_name]]), 
                 collapse = "\t")
})

# Write to a single GMT file
writeLines(unlist(gmt_lines), con = "./marker_genes/marker_sigs_filtered.gmt")