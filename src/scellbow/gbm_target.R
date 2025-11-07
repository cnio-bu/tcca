library(Seurat)
library(BPCells)
library(clustree)
library(tidyverse)
library(dplyr)

setwd("/storage/scratch01/shared/projects/bc-meta/")

bc <- readRDS("beyondcell/beyondcell_pancancer.rds")

# Subset GBM samples from gbm_nourhan_abdelfattah study
bc_gbm <- subset(bc, subset = study == "gbm_nourhan_abdelfattah")

# Write h5ad object from BPCells bc matrix as input to SCellBow
write_matrix_anndata_hdf5(mat = bc_gbm[["RNA"]]$data, 
                          path = "single_cell/scellbow/beyondcell/gbm_bcs_target.h5ad")

# Save the metadata for those cells
write.table(bc_gbm@meta.data,
            file = "single_cell/scellbow/beyondcell/gbm_metadata_target.tsv",
            row.names = TRUE,
            sep = "\t")

# Subset expression data
seu <- readRDS("single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
malignant <- subset(seu, subset = malignancy == TRUE)
colnames(malignant) <- paste0("c", c(1:ncol(malignant)))
write_matrix_anndata_hdf5(mat = malignant[["RNA"]]$counts.gbm_nourhan_abdelfattah,
                          path = "single_cell/scellbow/expression/gbm_expr_target.h5ad")

# Save the metadata for those cells
write.table(malignant@meta.data[malignant$study == "gbm_nourhan_abdelfattah", ],
            file = "single_cell/scellbow/expression/gbm_metadata_target.tsv",
            row.names = TRUE,
            sep = "\t")


# SEURAT ANALYSIS OF GLIOBLASTOMA STUDY
setwd("/storage/scratch01/users/mgonzalezb/bc-meta/scellbow/")
counts <- malignant[["RNA"]]$counts.gbm_nourhan_abdelfattah
metadata <- malignant@meta.data[malignant$study == "gbm_nourhan_abdelfattah", ]

gbm <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata
)

saveRDS(seu, "gbm_nourhan_abdelfattah.rds")

gbm <- readRDS("gbm_nourhan_abdelfattah.rds")

## Normalization, FindVariableFeatures and scaling
gbm <- NormalizeData(
    gbm,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)
gbm <- FindVariableFeatures(gbm, selection.method = "vst", n_features = 2000)

# Plot standarized variance of all genes
hvgs <- HVFInfo(gbm)
hvgs_sorted <- hvgs[order(-hvgs$variance.standardized), ]
png("plots/allgenes_variance.png", width = 10, height = 10, units = "in", res = 300)
plot(hvgs_sorted$variance.standardized,
  main = "Variance Explained by Features",
  ylab = "Standardized Variance",
  xlab = "Feature Rank"
)
dev.off()

# Plot variable features
top10 <- head(VariableFeatures(gbm), 10)
plot1 <- VariableFeaturePlot(gbm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png("plots/hvg_variance.png", width = 10, height = 10, units = "in", res = 300)
plot1 + plot2
dev.off()

# Select as highly variable genes those with a standarized variance higher than 2
top_hvgs <- hvgs_sorted %>%
  filter(variance.standardized >= 2)
# (around 1000 highly variable genes)
gbm <- FindVariableFeatures(gbm, selection.method = "vst", n_features = 1000)
gbm <- ScaleData(gbm)


# Run dimensionality reduction and clustering
gbm <- RunPCA(gbm, npcs = 100)

print(paste("PCA done", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

png("plots/elbow_plot.png", width = 10, height = 10, units = "in", res = 300)
ElbowPlot(gbm, ndims = 100)
dev.off()

gbm <- FindNeighbors(gbm, dims = 1:30)
gbm <- FindClusters(gbm, resolution = seq(0.1, 1, 0.1))

gbm <- RunUMAP(gbm, dims = 1:50)
gbm$patient_sample <- paste0(gbm$patient, "_", gbm$sample)

# Plot patient umap
gbm$patient <- factor(gbm$patient, levels = unique(gbm$patient))
patient.plot <- DimPlot(
  gbm,
  reduction = "umap",
  group.by = "patient",
  label = TRUE,
  label.size = 6
) + th
ggsave(
  "plots/umap_patient.png",
  plot = patient.plot,
  width = 10,
  height = 10,
  dpi = 300
)

# Rename patient_sample column and plot it
prefix <- gsub("_GSM\\d+", "", unique(gbm$patient_sample))
new_names <- ave(unique(gbm$patient_sample),
  prefix,
  FUN = function(x) paste0(unique(gsub("_GSM\\d+", "", x)), "_", seq_along(x))
)
mapping <- setNames(new_names, unique(gbm$patient_sample))
sample_renamed <- mapping[gbm$patient_sample]
names(sample_renamed) <- NULL
gbm$sample_renamed <- factor(sample_renamed, levels = unique(sample_renamed))

sample.plot <- DimPlot(
  gbm,
  reduction = "umap",
  group.by = "sample_renamed"
) + th
ggsave(
  "plots/umap_sample.png",
  plot = sample.plot,
  width = 12,
  height = 10,
  dpi = 300
)

# Plot tumor type on the umap
gbm$tumor_type_renamed <- ifelse(
  gbm$tumor_subtype == "diffuse astrocitoma",
  "AST",
  gbm$tumor_type
)

tumor_type.plot <- DimPlot(
  gbm,
  reduction = "umap",
  group.by = "tumor_type_renamed",
  label = TRUE,
  label.size = 6
) + th
ggsave(
  "plots/umap_tumor_type.png",
  plot = tumor_type.plot,
  width = 10,
  height = 10,
  dpi = 300
)

print(paste("UMAP done", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# We need to perform integration since each cluster correspond to a sample
gbm[["RNA"]] <- split(gbm[["RNA"]], f = gbm$patient_sample)
gbm <- NormalizeData(gbm)
gbm <- FindVariableFeatures(gbm)
gbm <- ScaleData(gbm)
gbm <- RunPCA(gbm)


gbm <- IntegrateLayers(
  object = gbm,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.harmony"
)

print(paste("Integration done", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

gbm <- JoinLayers(gbm)
gbm <- FindNeighbors(gbm, reduction = "integrated.harmony", dims = 1:30)
gbm <- FindClusters(gbm, resolution = 0.09)

gbm <- RunUMAP(
  gbm,
  reduction = "integrated.harmony",
  dims = 1:50,
  reduction.name = "umap.harmony"
)

saveRDS(gbm, "gbm_nourhan_abdelfattah_integrated.rds")


# Plot patient umap
gbm$patient <- factor(gbm$patient, levels = unique(gbm$patient))
patient.plot <- DimPlot(
  gbm,
  reduction = "umap.harmony",
  group.by = "patient"
) + th
ggsave(
  "plots/umap_patient_integrated.png",
  plot = patient.plot,
  width = 10,
  height = 10,
  dpi = 300
)

# Plot patient sample
sample.plot <- DimPlot(
  gbm,
  reduction = "umap.harmony",
  group.by = "sample_renamed"
) + th
ggsave(
  "plots/umap_sample_integrated.png",
  plot = sample.plot,
  width = 12,
  height = 10,
  dpi = 300
)

# Plot tumor type on the umap
tumor_type.plot <- DimPlot(
  gbm,
  reduction = "umap.harmony",
  group.by = "tumor_type_renamed",
) + th
ggsave(
  "plots/umap_tumor_type_integrated.png",
  plot = tumor_type.plot,
  width = 10,
  height = 10,
  dpi = 300
)

# Plot clusters on the integrated umap
clusters.plot <- DimPlot(
  gbm,
  reduction = "integrated.harmony",
  group.by = "RNA_snn_res.0.09"
)
ggsave(
  "plots/clusters_umap.png",
  plot = clusters.plot,
  width = 10,
  height = 10,
  dpi = 300
)

# Add VAR score to the plot
setwd("/storage/scratch01/shared/projects/bc-meta/")
var_scores <- read.table("single_cell/varscore/var_scores.tsv")

gbm@meta.data <- merge(gbm@meta.data, var_scores, by = "row.names", all = FALSE)
rownames(gbm@meta.data) <- gbm@meta.data$Row.names
gbm@meta.data$Row.names <- NULL


# Set the theme for all ggplots.
th <- theme(
  plot.title = element_text(hjust = 0.5, size = 17),
  legend.text = element_text(size = 15),
  legend.title = element_text(size = 15),
  axis.text.x = element_text(size = 15),
  axis.text.y = element_text(size = 15),
  axis.title.x = element_text(size = 15),
  axis.title.y = element_text(size = 15)
)

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/scellbow/")
a_scaled <- FeaturePlot(
  gbm, 
  reduction = "umap.harmony", 
  features = "A_scaled", 
  pt.size = 1.5
  ) +
  scale_color_viridis_c(option = "viridis", name = "A score") +
  ggtitle("Aggressiveness score (A scores)") + th

ggsave(
  "plots/ascore_umap.png",
  plot = a_scaled,
  width = 10,
  height = 10,
  dpi = 300
)

vr_scaled <- FeaturePlot(
  gbm,
  reduction = "umap.harmony",
  features = "VR_scaled",
  pt.size = 1.5
) +
  scale_color_viridis_c(option = "viridis", name = "VR score") +
  ggtitle("Vulnerability/Resistance score (VR scores)") + th

ggsave(
  "plots/vrscore_umap.png",
  plot = vr_scaled,
  width = 10,
  height = 10,
  dpi = 300
)


varsubs_scaled <- FeaturePlot(
  gbm,
  reduction = "umap.harmony",
  features = "VARscore_scaled_subs",
  pt.size = 1.5
) +
  scale_color_viridis_c(option = "viridis", name = "VAR score") +
  ggtitle("VAR score (A-VR score)") + th

ggsave(
  "plots/varsubs_umap.png",
  plot = varsubs_scaled,
  width = 10,
  height = 10,
  dpi = 300
)

varsum_scaled <- FeaturePlot(
  gbm,
  reduction = "umap.harmony",
  features = "VARscore_scaled_sum",
  pt.size = 1.5
) +
  scale_color_viridis_c(option = "viridis", name = "VAR score") +
  ggtitle("VAR score (A+VR score)") + th

ggsave(
  "plots/varsum_umap.png",
  plot = varsum_scaled,
  width = 10,
  height = 10,
  dpi = 300
)

varsum_scaled <- FeaturePlot(
  gbm,
  reduction = "umap.harmony",
  features = "VARscore_scaled_sum",
  pt.size = 1.5
) +
  scale_color_viridis_c(option = "viridis", name = "VAR score") +
  ggtitle("VAR score (A+VR score)") + th

ggsave(
  "plots/varsum_umap.png",
  plot = varsum_scaled,
  width = 10,
  height = 10,
  dpi = 300
)


