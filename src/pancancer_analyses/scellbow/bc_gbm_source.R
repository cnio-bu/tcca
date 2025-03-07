library(Seurat)
library(dplyr)
library(tidyverse)
library(openxlsx)
library(SingleCellExperiment)
library(zellkonverter)

setwd("/home/lmgonzalezb/Documents/bc-meta/SCellBow/GBM/source/")

# The source study not included in TCCA is Wang L, et al (2019).Cancer Discov. 
# doi: 10.1158/2159-8290.CD-19-0329. GEO: GSE138794

# Read scRNA-seq data from samples
samples <- list.dirs(".")[2:length(list.dirs())]

gbm <- lapply(samples, function(sample){
  gbm_data <- Read10X(data.dir = sample, gene.column = 1)
  gbm <- CreateSeuratObject(counts = gbm_data, 
                             project = sample, 
                             min.cells = 3, 
                             min.features = 200)
  return(gbm)
})


gbm <- merge(gbm[[1]], gbm[2: length(gbm)], add.cell.ids = gsub("\\./", "", samples))


### Preprocessing
gbm <- PercentageFeatureSet(gbm, pattern = "^MT-", col.name = "percent.mt")
gbm <- PercentageFeatureSet(gbm, pattern = "^RP[SL]", col.name = "percent.ribo")

gbm_subset <- subset(x = gbm, subset = (percent.mt <= 10) &
                       (nFeature_RNA >= 500 & nFeature_RNA <= 7000) &
                       (nCount_RNA > 500) & (percent.ribo <= 40)
)

this_counts <- GetAssayData(gbm_subset, slot = "counts")
nonzero_genes <- this_counts > 0

# Keep genes whose expression is found in at least 2% of the cells in the study
sample_cell_cutoff <- round(ncol(gbm_subset) / 100 * 2, digits = 0)
genes_to_keep <- Matrix::rowSums(nonzero_genes) >= sample_cell_cutoff

gbm_subset <- CreateSeuratObject(
  counts = this_counts[genes_to_keep, ],
  meta.data = gbm_subset@meta.data
)

# Remove cells from the same sample with duplicated names
cell_names <- gsub("_[0-9]$", "", colnames(gbm_subset))
cell_names_dup <- cell_names[duplicated(cell_names)]

cells_keep <- grep(paste0(cell_names_dup, collapse = "|"), 
                   colnames(gbm_subset), 
                   value = TRUE, 
                   invert = TRUE)

gbm_subset <- subset(gbm_subset, cells = cells_keep)
gbm_subset <- RenameCells(gbm_subset, 
                          new.names = gsub("_[0-9]$", "", 
                                           colnames(gbm_subset)))

# Add the cell type annotation
cell_type <- data.table::fread("GSE138794_scRNA_Seq_cell_types.txt.gz", header = FALSE) %>%
  as.data.frame()
rownames(cell_type) <- cell_type$V1
cell_type$V1 <- NULL
colnames(cell_type) <- c("cell_type")

common_cells <- colnames(gbm_subset)[colnames(gbm_subset) %in% rownames(cell_type)]
gbm_subset <- subset(gbm_subset, cells = common_cells)

gbm_subset$cell_type <- cell_type[colnames(gbm_subset),]
gbm_subset$sample <- gsub("^\\./", "", gbm_subset$orig.ident)

### Perform normalization, scaling
gbm_subset <- NormalizeData(gbm_subset,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000
)
gbm_subset <- FindVariableFeatures(gbm_subset, selection.method = "vst")
gbm_subset <- ScaleData(gbm_subset, features = rownames(gbm_subset))


# Run dimensionality reduction and clustering
gbm_subset <- RunPCA(gbm_subset, 
                     features = VariableFeatures(object = gbm_subset), 
                     npcs = 50)
ElbowPlot(gbm_subset, ndims = 50)

gbm_subset <- FindNeighbors(gbm_subset, dims = 1:30)
gbm_subset <- FindClusters(gbm_subset, resolution = 0.5)

gbm_subset <- RunUMAP(gbm_subset, dims = 1:30)

DimPlot(gbm_subset, reduction = "umap")
DimPlot(gbm_subset, reduction = "umap", group.by = "sample")
DimPlot(gbm_subset, reduction = "umap", group.by = "cell_type")

# Subset only malignant cells
gbm_malignant <- subset(gbm_subset, subset = cell_type == "Neoplastic_cell")


## Beyondcell
library(beyondcell)
library("beyondcell")
library("Seurat")
set.seed(1)

# Generate geneset object with one of the ready to use signature collections.
gs <- GetCollection(SSc, include.pathways = FALSE)
bc <- bcScore(gbm_malignant, gs, expr.thres = 0.1) 

bc@normalized[is.na(bc@normalized)] <- 0
bc <- bcRecompute(bc, slot = "normalized")
bc <- bcRegressOut(bc = bc, vars.to.regress = c("nFeature_RNA"))

# Save the Beyondcell object
saveRDS(bc, "bc_gbm_malignant.rds")

# Save the normalized beyondcell scores as a table
write.table(bc@normalized, 
            file = "gbm_mtx_source.tsv",
            sep = "\t")

# Save the metadata for those cells
write.table(bc@meta.data,
            file = "gbm_metadata_source.tsv",
            row.names = TRUE,
            sep = "\t")