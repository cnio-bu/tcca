library(BPCells)
library(Seurat)
library(clustree)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(reshape2)
library(ggpubr)
library(ComplexHeatmap)
library(ggalluvial)

setwd(
  "/storage/scratch01/users/mgonzalezb/bc-meta/studywise_analyses/perception_newmodels/"
)
# Load full perception matrix for this study
pc <- open_matrix_dir(dir = "v5/aml_sander_lambo")
pc_drugs <- readRDS("drug_models/new_drug_models.rds")
rownames(pc) <- names(pc_drugs)

# Load metadata for this samples
setwd("/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell_vs_perception/aml/")
seu <- readRDS("integrated_aml_sander_lambo.rds")
metadata <- seu@meta.data

colnames(pc) <- paste0("c", c(1:ncol(pc)))
rownames(metadata) <- paste0("c", c(1:ncol(pc)))

# Select common drugs
# drugs_bc <- readRDS("/drug_models/drugs_bc.rds")
# common_drugs <- intersect(rownames(pc), tolower(unique(drugs_bc)))
# pc <- pc[common_drugs, ]

# Create seurat object with perception score for clustering
seu <- Seurat::CreateSeuratObject(counts = pc,
                                  assay = "perception",
                                  meta.data = metadata)

seu[["perception"]]$data <- seu[["perception"]]$counts

# Masking so that variable flags get set.
VariableFeatures(seu) <- rownames(seu)

seu <- ScaleData(seu,
                 do.scale = TRUE,
                 do.center = TRUE,
                 scale.max = 20)

seu <- RunPCA(seu, reduction.name = "pca", npcs = 50)

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
  filename = "perception_clustering/new_models/all_drugs/loadings_lose_heats_33.png",
  heats,
  width = 29,
  height = 16,
  dpi = 300
)

elbow <- ElbowPlot(seu, ndims = 50)
ggsave(
  filename = "perception_clustering/new_models/all_drugs/elbow_plot_seu.png",
  elbow,
  width = 16,
  height = 9,
  dpi = 100
)

# Compute clusters based on reduced dimensions
seu <- FindNeighbors(object = seu,
                     dims = 1:50,
                     reduction = "pca")

seu <- FindClusters(
  seu,
  resolution = seq(0.1, 1, 0.1),
  # method = "igraph",
  # algorithm = 2,
  random.seed = 120394,
  group.singletons = TRUE
)

seu <- RunUMAP(seu, dims = 1:50)

clustree_tc <- clustree(seu@meta.data[, grep("perception_snn_res", colnames(seu@meta.data))], prefix = "perception_snn_res.")
ggsave(
  "perception_clustering/new_models/all_drugs/clustree.png",
  clustree_tc,
  width = 8,
  height = 12
)

for (i in seq(0.1, 1, 0.1)) {
  dimplot_clusters <- DimPlot(seu, group.by = paste0("perception_snn_res.", i))
  ggsave(
    paste0(
      "perception_clustering/new_models/all_drugs/umap_clusters",
      i,
      ".png"
    ),
    dimplot_clusters,
    width = 10,
    height = 8
  )
}

seu$Sample_type <- factor(sub("^[^_]*_([^_]*)_.*$", "\\1", seu$sample),
                          levels = c("Dx", "Rem", "Rel"))

dimplot_clusters <- DimPlot(seu, group.by = "patient")
ggsave(
  "perception_clustering/new_models/all_drugs/umap_patients.png",
  dimplot_clusters,
  width = 10,
  height = 8
)

dimplot_clusters <- DimPlot(seu, group.by = "Sample_type")
ggsave(
  "perception_clustering/new_models/all_drugs/umap_sample.png",
  dimplot_clusters,
  width = 10,
  height = 8
)



## Plot heatmap with therapeutic clusters
seu <- FindVariableFeatures(object = seu,
                            selection.method = "vst",
                            nfeatures = 79)
seu <- Seurat::SketchData(
  object = seu,
  assay = "perception",
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch_5k"
)

saveRDS(seu,
        "perception_clustering/new_models/all_drugs/seu_pc_clustering.rds")
print("File saved!")

mat <- as.matrix(seu[["sketch_5k"]]$data)
scaled_mat <- t(scale(
  x = t(-mat),
  center = TRUE,
  scale = TRUE
))


source("~/bc-meta/figures/TCCA_palette.R")

# Add treatment information
clinical <- read.xlsx("suplemmentary_aml_sander_lambo.xlsx", sheet = "Single cell Cohort") %>%
  mutate(patient = Patient.identifier) %>%
  select(patient, Treatment) %>%
  distinct()

seu@meta.data <- seu@meta.data %>%
  left_join(clinical, by = "patient")
rownames(seu@meta.data) <- paste0("c", seq_along(colnames(seu)))

top_annot_df <- seu@meta.data[colnames(mat), c("Sample_type",
                                               "Treatment",
                                               "Treatment_Outcome",
                                               "perception_snn_res.0.1")]
colnames(top_annot_df)[4] <- "Perception_clusters"
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
  "Perception_clusters" = c(
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
      preferred.drug.names == "MEVASTATIN" ~ "Other",!(preferred.drug.names %in% toupper(rownames(mat))) ~ "Other",
      TRUE ~ collapsed.MoAs
    )
  ) %>%
  distinct() %>%
  as.data.frame()

drugs <- drugs[tolower(drugs$preferred.drug.names) %in% rownames(mat), ] %>%
  select(preferred.drug.names, collapsed.MoAs) %>%
  distinct()
rownames(drugs) <- drugs$preferred.drug.names

MoAs <- drugs[toupper(rownames(mat)), "collapsed.MoAs"]
MoAs[is.na(MoAs)] <- "Other"
MoAs <- data.frame("Mechanism.of.action" = MoAs)
rownames(MoAs) <- rownames(mat)

moa_pals <- list("Mechanism.of.action" = MoAs_colors)


right_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = MoAs,
  which = "row",
  col = moa_pals,
  show_annotation_name = FALSE
)

png(
  file = "perception_clustering/new_models/all_drugs/heatmap_perception_clusters.png",
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
  column_order = rownames(top_annot_df[order(top_annot_df$`Perception_clusters`), ]),
  cluster_columns = FALSE,
  cluster_column_slices = TRUE,
  show_column_dend = FALSE,
  column_split = top_annot_df$`Perception_clusters`,
  clustering_distance_columns = "pearson",
  clustering_distance_rows = "pearson",
  show_column_names = FALSE,
  row_labels = toupper(rownames(mat)),
  show_row_names = TRUE,
  column_names_rot = 45,
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_side = "top",
  column_title = NULL,
  heatmap_legend_param = list(title = "Viability score"),
  heatmap_width = unit(14, "in"),
  heatmap_height = unit(8, "in")
)

ht_opt(
  "ANNOTATION_LEGEND_PADDING" = unit(1, "cm"),
  "HEATMAP_LEGEND_PADDING" = unit(1, "cm"),
  "legend_gap" = unit(1, "cm")
)
draw(heat, annotation_legend_side = "top")
dev.off()



# # Compare Beyondcell clustering with PERCEPTION clustering
# library(ggalluvial)
# library(ggplot2)

# seu_pc <- readRDS("perception_clustering/seu_pc_clustering.rds")
# seu_bc <- readRDS("beyondcell_clustering/seu_bc_clustering.rds")

# df <- data.frame("beyondcell_clusters" = seu_bc$beyondcell_clusters, "perception_clusters" = seu_pc$perception_clusters)
# write.table(df, "data_sankey.tsv")

# sankey <- ggplot(df, aes(axis1 = beyondcell_clusters, axis2 = perception_clusters)) +
#     geom_alluvium(aes(fill = beyondcell_clusters)) +
#     geom_stratum() +
#     geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#     theme_minimal() +
#     labs(title = "Sankey Plot for Cluster Resolutions",
#         x = "Tool",
#         y = "Count")
# ggsave(
#   "sankey_beyondcell_perception_clusters.png",
#   sankey,
#   dpi = 300,
#   height = 8,
#   width = 8
# )