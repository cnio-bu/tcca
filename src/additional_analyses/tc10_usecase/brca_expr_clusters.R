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
brca <- FindClusters(brca, resolution = seq(0.1, 1, 0.1))

brca <- RunUMAP(brca, dims = 1:40)
DimPlot(brca, reduction = "umap", group.by = "sample") + NoLegend()


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
clustree <- clustree(brca@meta.data[, grep("integrated_snn_res.", 
                                          colnames(brca@meta.data))], 
                     prefix = "integrated_snn_res.")
ggsave(
  "plots/clustree_integration.png",
  plot = clustree,
  dpi = 300,
  height = 7,
  width = 7
)

brca <- RunUMAP(
  brca,
  reduction = "integrated.harmony",
  dims = 1:40,
  reduction.name = "umap.harmony"
)
saveRDS(brca, "seu_brca_harmony.rds")

# To then use beyondcell for subclone level expression data, we need to transform
# the object to Seurat version 4
seurat_v4 <- JoinLayers(brca)

# Guardar
saveRDS(seurat_v4, "seu_brca_harmony_v4.rds")