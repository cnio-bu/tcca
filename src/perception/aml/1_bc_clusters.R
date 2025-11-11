library(BPCells)
library(Seurat)
library(clustree)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(reshape2)
library(ggpubr)
library(ComplexHeatmap)

setwd("/storage/scratch01/shared/projects/bc-meta")

# Load full perception matrix for this study
bc <- open_matrix_dir(dir = "./beyondcell/studywise_bpcells/aml_sander_lambo")
bc <- as.matrix(bc)
bc[is.na(bc)] <- 0
bc <- as(bc, "sparseMatrix")

# Load metadata for this samples
setwd("/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell_vs_perception/aml/")
seu <- readRDS("integrated_aml_sander_lambo.rds")
metadata <- seu@meta.data

colnames(bc) <- paste0("c", c(1:ncol(bc)))
rownames(metadata) <- paste0("c", c(1:ncol(bc)))

# Leave drugs common to perception and cells analyzed in previous pancancer beyondcell analysis
# pc <- open_matrix_dir("/storage/scratch01/users/mgonzalezb/bc-meta/studywise_analyses/perception_newmodels/full_mat_perception")

# load("drugInfo.RData")
# common_drugs <- drugInfo$IDs %>%
#     filter(collections == "SSc" & preferred.drug.names %in% toupper(rownames(pc))) %>%
#     pull(IDs)

# bc <- bc[common_drugs, ]


# Create seurat object with perception score for clustering
seu <- Seurat::CreateSeuratObject(counts = bc,
                                  assay = "beyondcell",
                                  meta.data = metadata)

seu[["beyondcell"]]$data <- seu[["beyondcell"]]$counts

# Masking so that variable flags get set.
VariableFeatures(seu) <- rownames(seu)

seu <- FindVariableFeatures(object = seu,
                            selection.method = "disp",
                            nfeatures = 300)

seu <- ScaleData(seu,
                 do.scale = TRUE,
                 do.center = TRUE,
                 scale.max = 20)

seu <- RunPCA(seu, reduction.name = "pca", npcs = 100)

## QC
heats <- DimHeatmap(
  seu,
  dims = 1:33,
  cells = 500,
  balanced = TRUE,
  slot = "scale.data",
  ncol = 3,
  nfeatures = 10,
  fast = FALSE,
  combine = TRUE
)

ggsave(
  filename = "beyondcell_clustering/all_drugs/loadings_lose_heats_33.png",
  heats,
  width = 29,
  height = 16,
  dpi = 300
)

elbow <- ElbowPlot(seu, ndims = 43)
ggsave(
  filename = "beyondcell_clustering/all_drugs/elbow_plot_seu.png",
  elbow,
  width = 16,
  height = 9,
  dpi = 100
)

# Compute clusters based on reduced dimensions
seu <- FindNeighbors(object = seu,
                     dims = 1:40,
                     reduction = "pca")

seu <- FindClusters(
  seu,
  resolution = seq(0.1, 1, 0.1),
  # method = "igraph",
  # algorithm = 2,
  random.seed = 120394,
  group.singletons = TRUE
)

seu <- RunUMAP(seu, dims = 1:40)

clustree_tc <- clustree(seu@meta.data[, grep("beyondcell_snn_res", colnames(seu@meta.data))], prefix = "beyondcell_snn_res.")
ggsave(
  "beyondcell_clustering/all_drugs/clustree.png",
  clustree_tc,
  width = 8,
  height = 12
)

for (i in seq(0.1, 1, 0.1)) {
  dimplot_clusters <- DimPlot(seu, group.by = paste0("beyondcell_snn_res.", i))
  ggsave(
    paste0("beyondcell_clustering/all_drugs/umap_clusters", i, ".png"),
    dimplot_clusters,
    width = 10,
    height = 8
  )
}


seu$Sample_type <- factor(sub("^[^_]*_([^_]*)_.*$", "\\1", seu$sample),
                          levels = c("Dx", "Rem", "Rel"))






## Plot heatmap with therapeutic clusters

seu <- Seurat::SketchData(
  object = seu,
  assay = "beyondcell",
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch_5k"
)

saveRDS(seu,
        "beyondcell_clustering/all_drugs/seu_bc_clustering.rds")

sketched_mat <- as.matrix(seu[["sketch_5k"]]$data)

scaled_mat <- t(scale(
  x = t(sketched_mat),
  center = TRUE,
  scale = TRUE
))
top_rv <- matrixStats::rowVars(scale(x = sketched_mat, center = TRUE, scale = TRUE))
top_rv <- top_rv[top_rv >= 2]

source("~/bc-meta/figures/TCCA_palette.R")

# Add treatment information
clinical <- read.xlsx("suplemmentary_aml_sander_lambo.xlsx", sheet = "Single cell Cohort") %>%
  mutate(patient = Patient.identifier) %>%
  select(patient, Treatment) %>%
  distinct()

seu@meta.data <- seu@meta.data %>%
  left_join(clinical, by = "patient")
rownames(seu@meta.data) <- paste0("c", seq_along(colnames(seu)))

top_annot_df <- seu@meta.data[colnames(scaled_mat), c("Sample_type",
                                                      "Treatment",
                                                      "Treatment_Outcome",
                                                      "beyondcell_snn_res.0.1")]
colnames(top_annot_df)[4] <- "Beyondcell_clusters"
pals <- list(
  "Sample_type" = c(
    "Rel" = "#C10044",
    "Dx" = "#F0BFD0",
    "Rem" = "#706695"
  ),
  "Treatment" = c("Arm B" = "#427394", "Arm C" = "#F78C1F"),
  "Treatment_Outcome" = c(
    "Relapsed" = "#6ED1BC",
    "Censored" = "#D18B6E"
  ),
  "Beyondcell_clusters" = c(
    "0"  = "#FFBD71",
    "1"  = "#FFA72C",
    "2"  = "#FE7B47",
    "3"  = "#D05B61",
    "4"  = "#FFB0EA",
    "5"  = "#B46F9C",
    "6"  = "#A52390",
    "7"  = "#6567BD",
    "8"  = "#406792",
    "9"  = "#369CBB",
    "10" = "#90E2ED",
    "11" = "#43978D",
    "12" = "#A7D676",
    "13" = "#D2E295",
    "14" = "#DBDF00",
    "15" = "#FFE364"
  )
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df =  top_annot_df,
  which = "column",
  col = pals,
  annotation_name_side = "left",
  annotation_name_rot = 0
)


## get drug names
drugs <- data.table::fread("~/bc-meta/reference/final_moas - Collapsed.tsv") %>%
  select(IDs, preferred.drug.names, collapsed.MoAs) %>%
  mutate(
    collapsed.MoAs = case_when(
      preferred.drug.names == "VANDETANIB" ~ "VEGFR inhibitor",
      preferred.drug.names == "DASATINIB" ~ "Kinase inhibitor",
      preferred.drug.names == "RIGOSERTIB" ~ "Other",
      preferred.drug.names == "TALAZOPARIB" ~ "Other",
      preferred.drug.names == "SORAFENIB" ~ "Kinase inhibitor",
      preferred.drug.names == "SARACATINIB" ~ "Other",
      preferred.drug.names == "PITAVASTATIN" ~ "Other",
      preferred.drug.names == "MEVASTATIN" ~ "Other",
      TRUE ~ collapsed.MoAs
    )
  ) %>%
  distinct() %>%
  as.data.frame()

drugs <- drugs[drugs$IDs %in% rownames(scaled_mat), ]
rownames(drugs) <- drugs$IDs

MoAs <- drugs[rownames(scaled_mat), c("IDs", "collapsed.MoAs")]
MoAs <- as.data.frame(MoAs$collapsed.MoAs)
colnames(MoAs) <- "Mechanism of action"

moa_pals <- list("Mechanism of action" = MoAs_colors)


right_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = MoAs,
  which = "row",
  col = moa_pals,
  show_annotation_name = FALSE
)

png(
  file = "beyondcell_clustering/all_drugs/heatmap_beyondcell_with_tcs.png",
  res = 500,
  width = 19,
  height = 14,
  units = "in"
)

# Customize legends for tumor site
heat <- ComplexHeatmap::Heatmap(
  mat = scaled_mat,
  #mat = t(sketched_mat),
  right_annotation = right_annotation,
  top_annotation = top_annotation,
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  row_split = 5,
  row_title = NULL,
  column_order = rownames(top_annot_df[order(top_annot_df$`Beyondcell_clusters`), ]),
  cluster_columns = FALSE,
  cluster_column_slices = TRUE,
  show_column_dend = FALSE,
  column_split = top_annot_df$`Beyondcell_clusters`,
  clustering_distance_columns = "pearson",
  clustering_distance_rows = "pearson",
  show_column_names = FALSE,
  row_labels = paste0(drugs[rownames(scaled_mat), "preferred.drug.names"]),
  show_row_names = TRUE,
  column_names_rot = 45,
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_side = "top",
  column_title = NULL,
  heatmap_legend_param = list(title = "BCS score"),
  heatmap_width = unit(14, "in"),
  heatmap_height = unit(10, "in")
)

ht_opt(
  "ANNOTATION_LEGEND_PADDING" = unit(1, "cm"),
  "HEATMAP_LEGEND_PADDING" = unit(1, "cm"),
  "legend_gap" = unit(1, "cm")
)
draw(heat, annotation_legend_side = "top")
dev.off()

bc <- readRDS("beyondcell_clustering/all_drugs/seu_bc_clustering.rds")
pc <- readRDS("aml/perception_clustering/new_models/all_drugs/seu_pc_clustering.rds")

df <- data.frame(
  "beyondcell_clusters" = bc$beyondcell_snn_res0.1,
  "perception_clusters" = pc$perception_snn_res.0.05
)
write.table(df, "aml_clusters_all_drugs.tsv")

# Save objects with common drugs
bc <- readRDS("beyondcell_clustering/common_drugs/seu_bc_clustering.rds")
pc <- readRDS("perception_clustering/new_models/common_drugs/seu_pc_clustering.rds")
df <- data.frame(
  "beyondcell_clusters" = bc$beyondcell_snn_res.0.1,
  "perception_clusters" = pc$RNA_snn_res.0.1
)
write.table(df, "aml_clusters_common_drugs.tsv")