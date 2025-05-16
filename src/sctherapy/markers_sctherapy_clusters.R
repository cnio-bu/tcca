library(BPCells)
library(Seurat)
library(dplyr)
library(tidyverse)
library(presto)

setwd("/storage/scratch01/shared/projects/bc-meta/")
seu <- readRDS("./single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")

# Leave only malignant cells with clonal information from SCEVAN
full_metadata <- read.table("./functional_nmf/subclone_wise/metadata_subclone_annot.tsv")

# Subset the malignant cells with clonal information
seu <- subset(seu, cells = rownames(full_metadata))
seu <- AddMetaData(seu, metadata = full_metadata)

# Add the scTherapy clusters to the metadata
clusters <- readRDS("single_cell/sctherapy/results/speclustering_reordered.rds")
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

saveRDS(seu, "single_cell/sctherapy/results/seu_before_cluster_markers.rds")
# Find differentially expressed features across scTherapy clusters
seu <- JoinLayers(seu)
DefaultAssay(seu) <- "RNA"
Idents(seu) <- "clusters"
seu.markers <- FindAllMarkers(seu, only.pos = TRUE)

# Save top 100 genes per cluster as a table
seu.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 2) %>%
    ungroup() -> all_markers

all_markers <- all_markers %>%
    select(cluster, gene) %>%
    group_by(cluster) %>%
    mutate(row = row_number()) %>%
    ungroup() %>%
    pivot_wider(names_from = cluster, values_from = gene) %>%
    select(-row) %>%
    as.data.frame()

write.table(as.data.frame(seu.markers), "single_cell/sctherapy/results/marker_genes/table_markers_stats.tsv", 
            sep = "\t", row.names = FALSE)
all_markers <- all_markers[as.character(1:10)]
write.table(all_markers, "single_cell/sctherapy/results/marker_genes/all_markers.tsv", 
            sep = "\t", row.names = FALSE)


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
    dplyr::filter(avg_log2FC > 2) %>%
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


png("./single_cell/sctherapy/results/marker_genes/heatmap_top10_cluster.png", 
    width = 12, 
    height = 10, 
    units = "in",
    res = 500)
DoHeatmap(seu, features = top10$gene, group.by = "clusters") + NoLegend()
dev.off()

