library(tidyverse)
library(BPCells)
library(Seurat)
library(clustree)


## Set options
options(future.globals.maxSize = 20 * 1024 ^ 3)
options(Seurat.object.assay.version = 'v5')

setwd("/storage/scratch01/users/mgonzalezb/bc-meta/studywise_analyses/perception")

# #Load full perception mat and drug models
mat <- open_matrix_dir(dir = "./full_mat_perception")
drug_models <- readRDS("./drug_models/FDA_approved_drugs_models.rds")
rownames(mat) <- names(drug_models)

# Load full database
meta.data_full_clinical <- read_tsv(file = "metadata_with_clinical.tsv")  %>% as.data.frame()

meta.data_full_clinical$new_cell_id <-  paste0("c", c(1:nrow(meta.data_full_clinical)))
rownames(meta.data_full_clinical) <- meta.data_full_clinical$new_cell_id

colnames(mat) <- paste0("c", c(1:ncol(mat)))

meta.data_full_clinical$new_cell_id <-  paste0("c", c(1:nrow(meta.data_full_clinical)))
rownames(meta.data_full_clinical) <- meta.data_full_clinical$new_cell_id

seu <- Seurat::CreateSeuratObject(
  counts = mat,
  assay = "RNA",
  project = "perception_pancancer",
  meta.data = meta.data_full_clinical
)

seu[["RNA"]]$data <- seu[["RNA"]]$counts


# Load the beyondcell pancancer object to select only cells where BCS was
# computed in order to compare with perception results
bc <- readRDS(
  "/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell/results/beyondcell_pancancer_final_res.Rds"
)
seu <- subset(seu, cells = colnames(bc))

## Use the full dataset as variable feat. We are talking 581. Let PCA summarise.
seu <- FindVariableFeatures(seu, selection.method = "vst")

## Masking so that variable flags get set.
VariableFeatures(seu) <- rownames(seu)

## Go for sketch, it's faster
seu <- SketchData(
  object = seu,
  ncells = 50000,
  method = "LeverageScore",
  over.write = TRUE,
  sketched.assay = "sketch_50k",
  seed = 120394
)

DefaultAssay(seu) <- "sketch_50k"

seu <- FindVariableFeatures(object = seu,
                            selection.method = "disp",
                            nfeatures = 300)

seu <- ScaleData(seu,
                 do.scale = TRUE,
                 do.center = TRUE,
                 scale.max = 20)

seu <- RunPCA(seu, reduction.name = "pca", npcs = 100)

## QC
vimdizload <- VizDimLoadings(seu, dims = 1:2, reduction = "pca")

ggsave(
  filename = "results/plots_qc/vim_loadings_pca.png",
  plot = vimdizload,
  height = 9,
  width = 16,
  dpi = 100
)

seu_pca <- DimPlot(seu, reduction = "pca")

ggsave(
  filename = "results/plots_qc/seu_pca_sketched_50.png",
  plot = seu_pca,
  height = 9,
  width = 16,
  dpi = 100
)

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
  filename = "results/plots_qc/loadings_lose_heats_33.png",
  heats,
  width = 29,
  height = 16,
  dpi = 300
)

heats <- DimHeatmap(
  seu,
  dims = 1:10,
  cells = 500,
  balanced = TRUE,
  slot = "scale.data",
  ncol = 2,
  nfeatures = 10,
  fast = FALSE,
  combine = TRUE
)

ggsave(
  filename = "results/plots_qc/loadings_lose_heats_10.png",
  heats,
  width = 16,
  height = 9,
  dpi = 300
)

elbow <- ElbowPlot(seu, ndims = 43)
ggsave(
  filename = "results/plots_qc/elbow_plot_seu.png",
  elbow,
  width = 16,
  height = 9,
  dpi = 100
)


seu <- FindNeighbors(
  object = seu,
  dims = 1:40,
  reduction = "pca",
  k.param = 300
)
seu <- FindClusters(
  seu,
  resolution = seq(0.1, 1, 0.1),
  # method = "igraph",
  # algorithm = 2,
  random.seed = 120394,
  group.singletons = TRUE
  #cluster.name = "therapeutic_clusters_sketched_k.300_res.0.5"
)

seu <- RunUMAP(
  seu,
  dims = 1:40,
  return.model = T,
  umap.method = "uwot",
  #n.neighbors = 50,
  #metric = "correlation",
  #local.connectivity = 4
)

clusters_project <- setNames(
  paste0("sketch_50k_snn_res.", seq(0.1, 1, 0.1)),
  # Values
  paste0("therapeutic_clusters_k.300.res.", seq(0.1, 1, 0.1))  # Names
)

seu <- Seurat::ProjectData(
  object = seu,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch_50k",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:40,
  refdata = clusters_project,
)

write_tsv(seu@meta.data, "results/tcs.tsv")

seu <- Seurat::SketchData(
  object = seu,
  assay = "RNA",
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch_5k"
)

saveRDS(object = seu, file = "./results/perception_pancancer_FDA_approved.Rds")

## export  sketches
sketched_mat <- seu[["sketch_50k"]]$data
sketched_mat_5k <- seu[["sketch_5k"]]$data
write_matrix_dir(mat = sketched_mat,
                 dir = "./results/sketch_mat_perception",
                 overwrite = TRUE)
write_matrix_dir(mat = sketched_mat_5k,
                 dir = "./results/sketch_mat_perception_5k",
                 overwrite = TRUE)
