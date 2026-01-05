library(Seurat)
library(BPCells)
library(dplyr)
library(tidyverse)
library(clustree)
setwd("/storage/scratch01/shared/projects/bc-meta/beyondcell_immuno/brca_usecase")

# Load expression matrix of cells from BRCA patients
expr.mat <- readRDS("input_data/expr_sparsemat.rds")

# Filter metadata for BRCA patients
metadata <- read.table(
  "../../single_cell/seurat/v5/tcca_metadata.tsv",
  sep = "\t", header = TRUE
)

metadata_brca <- metadata %>%
  filter(malignancy == "True") %>%
  mutate(cell = row_number()) %>%
  filter(refined_tumor_type == "BRCA" & patient != "ccl")
rownames(metadata_brca) <- metadata_brca$cell

# Create Seurat object
brca <- CreateSeuratObject(counts = expr.mat, meta.data = metadata_brca)

# Perform Seurat analysis.
brca <- NormalizeData(brca)

brca <- FindVariableFeatures(brca)
top10 <- head(VariableFeatures(brca), 10)
plot1 <- VariableFeaturePlot(brca)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

brca <- ScaleData(brca)
brca <- RunPCA(brca, npcs = 100)
ElbowPlot(brca, ndims = 100)


brca <- FindNeighbors(brca, dims = 1:40)
brca <- FindClusters(brca, resolution = 0.5, cluster.name = "unintegrated_umap")

brca <- RunUMAP(brca, dims = 1:40)
cells_tcs <- rownames(brca@meta.data %>% filter(!is.na(scTherapy_cluster)))
source("/home/mgonzalezb/bc-meta/figures/TCCA_palette.R")
png("plots/dimplot_tcs_unintegrated.png", width = 5, height = 4, units = "in", res = 300)
DimPlot(
  brca, 
  reduction = "umap", 
  cells = cells_tcs, 
  group.by = "scTherapy_cluster", 
  cols = sctherapy_colors
  ) + 
  NoLegend()
dev.off()

# Large batch effect caused by different samples in the visualization
brca$study_sample <- paste0(brca$study, "_", brca$sample)
brca[["RNA"]] <- split(brca[["RNA"]], f = brca$study_sample)

# Preprocessing steps for each layer (sample) 
brca <- NormalizeData(brca)
brca <- FindVariableFeatures(brca)
brca <- ScaleData(brca)
brca <- RunPCA(brca)

brca <- IntegrateLayers(
  object = brca,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony"
)

brca <- JoinLayers(brca)
brca <- FindNeighbors(brca, reduction = "integrated.harmony", dims = 1:40)
brca <- FindClusters(brca, resolution = seq(0.1, 1, 0.1))

clustree(brca@meta.data[, grep("RNA_snn_res.", colnames(brca@meta.data))], 
         prefix = "RNA_snn_res.")

brca <- RunUMAP(
  brca,
  reduction = "integrated.harmony",
  dims = 1:40,
  reduction.name = "umap.harmony"
)
saveRDS(brca, "seu_brca_harmony.rds")


# Move to local directories
setwd("/Users/mariagb/Documents/bc_meta/brca_usecase")
brca <- readRDS("seu_brca_harmony.rds")
png("plots/dimplot_tcs_integrated.png", width = 5, height = 4, units = "in", res = 300)
DimPlot(
  brca, 
  reduction = "umap.harmony", 
  cells = cells_tcs, 
  group.by = "scTherapy_cluster",
  cols = sctherapy_colors
  )
dev.off()
