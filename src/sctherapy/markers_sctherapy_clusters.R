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
metadata_tcca <- read.table("./single_cell/seurat/tcca/tcca_metadata.tsv", header = TRUE)

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

saveRDS(seu, "seu_before_cluster_markers.rds")

# Find differentially expressed features across scTherapy clusters
seu <- JoinLayers(seu)
DefaultAssay(seu) <- "RNA"
Idents(seu) <- "clusters"
seu.markers <- FindAllMarkers(seu, only.pos = TRUE)
seu.markers.auc <- FindAllMarkers(seu, test.use = "roc", return.thresh = 0.6, only.pos = TRUE)
write.table(seu.markers.auc, "gene_markers_auc2.tsv", row.names = FALSE, col.names = TRUE)

# Plot proportion of samples per cluster
seu@meta.data <- seu@meta.data %>%
  left_join(select(metadata_tcca, cell, refined_tumor_type), by = "cell")

barplot <- ggplot(seu@meta.data, aes(x = clusters, fill = refined_tumor_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = tumor_type_colors) +
  labs(x = "Cancer type", y = "Cell fraction", fill = "Tumor type") +
  ggtitle("Tumor types of cells across clusters") +
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 14, margin = margin(t = 6)),
        axis.title.y = element_text(size = 14, margin = margin(r = 6)),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
ggsave("marker_genes/barplot_cell_cancertype.pdf", plot = barplot, height = 6, width = 10)


# Save top genes per cluster as a table
seu.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC >= log2(1.5) & pct.1 >= 0.5 & (pct.1 - pct.2) >= 0.20) %>%
  arrange(desc(avg_log2FC), desc(pct.1)) %>%
  ungroup() -> all_markers

all_markers <- all_markers %>%
  select(cluster, gene) %>%
  group_by(cluster) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = cluster, values_from = gene) %>%
  select(-row) %>%
  as.data.frame()

write.table(as.data.frame(seu.markers), "marker_genes/table_markers_stats.tsv", 
            sep = "\t", row.names = FALSE)
all_markers <- all_markers[as.character(1:10)]
write.table(all_markers, "marker_genes/all_markers.tsv", 
            sep = "\t", row.names = FALSE)


# Plot top markers in a heatmap
seu.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC >= log2(1.5) & pct.1 > 0.5 & (pct.1 - pct.2) >= 0.20) %>%
  arrange(desc(avg_log2FC), desc(pct.1)) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15

expr_mat <- seu[["sketch"]]$scale.data[top10$gene,]

metadata <- seu@meta.data[colnames(expr_mat), ]
metadata <- metadata %>%
  left_join(select(metadata_tcca, cell, refined_tumor_type), by = "cell") %>%
  select(cell, clusters, refined_tumor_type) %>%
  column_to_rownames(var = "cell")

metadata$clusters <- factor(metadata$clusters)
metadata$refined_tumor_type <- factor(metadata$refined_tumor_type)

pals <- list(
  "Clusters" = sctherapy_colors,
  "Tumor type" = tumor_type_colors
)
colnames(metadata) <- c("Clusters", "Tumor type")
col_anno <- ComplexHeatmap::HeatmapAnnotation(
  df = metadata,
  col = pals,
  annotation_name_side = "left",
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
  annotation_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "bold"),
                                 labels_gp = gpar(fontsize = 12),
                                 title_gap = unit(10, "mm")),
  show_legend = c(
    "Tumor type" = FALSE
  )
)

tumor_type_legend <- Legend(
  at = names(pals$`Tumor type`),
  legend_gp = gpar(fill = pals$`Tumor type`),
  title_gp = gpar(fontsize = 12, fontface = "bold"),
  labels_gp = gpar(fontsize = 12),
  ncol = 4,  # Split Group legend into 2 columns
  gap = unit(5, "mm"),
  title = "Tumor type"
)

heatmap_top10 <- Heatmap(
  as.matrix(expr_mat),
  name = "Expression",
  top_annotation = col_anno,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  column_split = metadata$Clusters,
  heatmap_width = unit(10, "in"),
  heatmap_height = unit(12, "in")
)

png(
  file = "marker_genes/heatmap_markers.png",
  res = 500,
  width = 14,
  height = 16,
  units = "in"
)
draw(heatmap_top10 , 
     annotation_legend_side = "top", 
     heatmap_legend_side = "right", 
     annotation_legend_list = list(tumor_type_legend))
dev.off()


# Run wilcoxauc to get the AUC and the wilcoxon statistics
expr_mat <- as(seu[["RNA"]]$data, "dgCMatrix")
markers.auc <- presto::wilcoxauc(expr_mat, seu$clusters)

write.table(markers, "gene_markers_auc.tsv", row.names = FALSE, col.names = TRUE)

markers.auc <- FindAllMarkers(seu, only.pos = TRUE)
markers.auc.top10 <- markers.auc %>%
  group_by(group) %>%
  dplyr::filter(logFC >= log(1.5) & pct_in > 50 & (pct_in - pct_out) >= 20) %>%
  arrange(desc(auc), desc(logFC)) %>%
  slice_head(n = 15) %>%
  ungroup()



# Plot top markers per cluster as a dotplot
DefaultAssay(seu) <- "RNA"
dot_markers <- DotPlot(
  seu, 
  features = unique(top15$gene),
  group.by = "clusters"
) +
  coord_flip()
ggsave("dot_markers.png", plot = dot_markers, height = 10, width = 10)

# Extract average expression and percentage expressed from the dot plot
df <- dot_markers$data
exp_mat <- df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

head(exp_mat)

# Extract the percentage of cells express a gene
percent_mat <-df %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame()

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

head(percent_mat)
range(percent_mat)
dim(exp_mat)
dim(percent_mat)

# Create color palette
col_fun <- circlize::colorRamp2(c(-2, 0, 2), magma(20)[c(20, 10, 1)])
max_r <- unit(0.4, "snpc")
cell_fun <- function(j, i, x, y, w, h, fill) {
  grid.rect(
    x = x, y = y, width = w, height = h,
    gp = gpar(col = NA, fill = NA)
  )
  grid.circle(
    x = x, y = y, r = percent_mat[i, j] / 100 * min(unit.c(w, h)),
    gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA)
  )
}

dot_legend <- Legend(
  title = "% Expressed",
  at = c(25, 50, 75, 100),
  type = "points",
  legend_gp = gpar(fill = "black"),
  pch = 21,
  size = unit(c(0.25, 0.50, 0.75, 1), "cm") # Adjust size to match circle scaling
)

# Add cluster annotation
column_ha <- HeatmapAnnotation(
  cluster_anno = as.factor(c(1:10)),
  col = list(cluster_anno = sctherapy_colors),
  na_col = "grey"
)

heatmap <- Heatmap(exp_mat,
                   heatmap_legend_param = list(title = "expression"),
                   column_title = "clustered dotplot",
                   col = col_fun,
                   rect_gp = gpar(type = "none"),
                   cell_fun = cell_fun,
                   row_names_gp = gpar(fontsize = 5),
                   border = "black",
                   top_annotation = column_ha
)
png("heatmap_test.png", width = 6, height = 8, res = 300, unit = "in")
draw(heatmap, annotation_legend_side = "right", annotation_legend_list = list(dot_legend))
dev.off()

# Add annotation to the heatmap
colnames(exp_mat)
library(RColorBrewer)
cluster_anno <- c("CD4T", "B", "CD4T", "Mono", "NK", "CD8T", "CD14_Mono", "DC", "Platelet")

column_ha <- HeatmapAnnotation(
  cluster_anno = as.factor(c(1:10)),
  col = list(cluster_anno = sctherapy_colors),
  na_col = "grey"
)

Heatmap(exp_mat,
        heatmap_legend_param = list(title = "expression"),
        column_title = "clustered dotplot",
        col = col_fun,
        rect_gp = gpar(type = "none"),
        cell_fun = cell_fun,
        row_names_gp = gpar(fontsize = 5),
        row_km = 4,
        border = "black",
        top_annotation = column_ha
)


# Select top genes and plot
seu.markers <- seu.markers %>%
  mutate(cluster_gene = paste0(cluster, ".", gene))
seu.markers.auc <- seu.markers.auc %>%
  mutate(cluster_gene = paste0(cluster, ".", gene))

top_genes <- seu.markers %>%
  left_join(select(seu.markers.auc, myAUC, avg_diff, power, cluster_gene), by = "cluster_gene") %>%
  group_by(cluster) %>%
  # dplyr::filter(avg_log2FC >= log2(1.5) & pct.1 > 0.5 & (pct.1 - pct.2) >= 0.20 & myAUC >= 0.60) %>%
  arrange(cluster, desc(myAUC), desc(avg_diff), desc(pct.1), desc(avg_log2FC)) %>%
  ungroup()

signatures <- top_genes %>%
  select(cluster, gene) %>%
  mutate(set_name = paste0("Cluster", stringr::str_pad(cluster, width = 2, pad = "0"))) %>%
  group_by(set_name) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  deframe()

names(signatures) <- paste0(names(signatures), "_UP")

# Convert to GMT lines
gmt_lines <- lapply(names(signatures), function(sig_name){
  sig <- paste(c(sig_name, "marker_genes_sctherapy_clusters", signatures[[sig_name]]), 
               collapse = "\t")
})

# Write to a single GMT file
writeLines(unlist(gmt_lines), 
           con = "./markers/marker_sigs_clusters_clean.gmt")


# png("heatmap_test.png", width = 7, height = 11, res = 300, unit = "in")
scCustomize::Clustered_DotPlot(
  seu,
  features = top15$gene, 
  group.by = "clusters",
  colors_use_idents = sctherapy_colors,
  legend_label_size = 8,
  legend_title_size = 8
)
dev.off()

genes_band <- read.table("marker_genes/genes_in_17q21.2.tsv", sep = "\t", header = TRUE)
c4 <- top_genes[top_genes$cluster == 4, "gene"]
common_c4 <- intersect(c4, genes_band$external_gene_name)














# Create signatures of genes with top25, top50, top75, top100, top125 and top150
get_top_genes <- function(df, n){
  df %>%
    filter(avg_log2FC > 2) %>%
    group_by(cluster) %>%
    slice_head(n = n) %>%
    ungroup() %>%
    select(cluster, gene) %>%
    mutate(set_name = paste0("Cluster", stringr::str_pad(cluster, width = 2, pad = "0"), "_top", n)) %>%
    group_by(set_name) %>%
    summarise(genes = list(unique(gene)), .groups = "drop") %>%
    deframe()
}


top25 <- get_top_genes(seu.markers, 25)
top50 <- get_top_genes(seu.markers, 50)
top75 <- get_top_genes(seu.markers, 75)
top100 <- get_top_genes(seu.markers, 100)
top125 <- get_top_genes(seu.markers, 125)
top150 <- get_top_genes(seu.markers, 150)

# Combine all
all_top <- c(top25, top50, top75, top100, top125, top150)
names(all_top) <- paste0(names(all_top), "_UP")

# Convert to GMT lines
gmt_lines <- lapply(names(all_top), function(sig_name){
  sig <- paste(c(sig_name, "marker_genes_sctherapy_clusters", all_top[[sig_name]]), 
               collapse = "\t")
})

# Write to a single GMT file
writeLines(unlist(gmt_lines), 
           con = "./single_cell/sctherapy/results/marker_genes/marker_sigs_clusters.gmt")

write.gmt <- function(df, file){
  
}
saveRDS(seu.markers, "single_cell/sctherapy/results/table_markers.rds")
saveRDS(seu, "single_cell/sctherapy/results/seu_cluster_markers.rds")

# Compute sketch of 5k cells
features_to_include <- unique(unlist(top150))
seu <- SketchData(
  object = seu,
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch",
  features = features_to_include
)

saveRDS(seu, "single_cell/sctherapy/results/seu_cluster_markers.rds")

# Create a heatmap with top marker genes
seu <-  ScaleData(seu, features = rownames(seu))

seu.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC >= log2(1.5) & pct.1 >= 0.5 & (pct.1 - pct.2) >= 0.20) %>%
  arrange(desc(avg_log2FC), desc(pct.1)) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# Add broad cancer types
cancer_type <- list(
  "Brain Cancer" = c("GBM", "MB", "OGD"),
  "Neuroblastic Tumors" = c("GNB", "NB"),
  "Blood Cancer" = c("ALL", "LAML", "CLL", "MM"),
  "Skin Cancer" = c("BCC", "SKCM", "SKSC", "SKAM", "UVM"),
  "Sarcoma/Soft Tissue Cancer" = c("SARC", "GIST", "MESO"),
  "Breast Cancer" = c("BRCA"),
  "Lung Cancer" = c("SCLC", "NSCLC", "LUAD", "LUSC", "LCLC", "PLEU"),
  "Ovarian Cancer" = c("OV"),
  "Colon/Colorectal Cancer" = c("COAD", "READ"),
  "Endometrial/Uterine Cancer" = c("CESC", "UCEC", "UCS"),
  "Liver/Biliary Cancer" = c("LIHC", "CHOL"),
  "Bladder Cancer" = c("BLCA"),
  "Head and Neck Cancer" = c("HNSC"),
  "Prostate Cancer" = c("PRAD"),
  "Kidney Cancer" = c("KRCC", "KTCC", "KIRC", "KIRCH"),
  "Esophageal Cancer" = c("ESCA", "ESCC"),
  "Pancreatic Cancer" = c("PAAD"),
  "Thyroid Cancer" = c("THCA"),
  "Gastric Cancer" = c("STAD"),
  "Miscellaneous Cancer" = c("MISC")
)

cancer_type <- enframe(cancer_type, name = "broad_cancer_type", value = "tumor_type") %>%
  unnest()

seu@meta.data <-  seu@meta.data %>%
  rownames_to_column(var = "cell_id") %>%
  left_join(cancer_type, by = "tumor_type") %>%
  column_to_rownames(var = "cell_id")


png("./single_cell/sctherapy/results/marker_genes/heatmap_top10_cancer_type.png", 
    width = 12, 
    height = 10, 
    units = "in",
    res = 500)
DoHeatmap(seu, features = top10$gene, group.by = "broad_cancer_type.x") + NoLegend()
dev.off()


# New signature cluster 10
markers_clust10 <- gene.markers %>%
  filter(cluster == 10 & avg_log2FC >= 0.6 & pct.1 >= 0.3 & pct.2 <= 0.1)

