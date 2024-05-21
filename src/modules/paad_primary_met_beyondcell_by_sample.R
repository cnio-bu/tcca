library(beyondcell)

mat <- readRDS("results/paad/normalized_seu_malignant_mat.rds")
meta.data <- read.table("results/paad/seu_metadata_malignants.tsv")

cells_to_keep <- meta.data[meta.data$patient == "P02", ]
mat <- mat[, rownames(cells_to_keep)]

#gs <- GetCollection(SSc, include.pathways = FALSE)
gs <- GenerateGenesets(
    x = "reference/drug_signatures_fold.gmt",
    perform.reversal = FALSE
)

bc <- bcScore(sc = mat, gs = gs, expr.thres = 0.1)
bc@normalized[is.na(bc@normalized)] <- 0
bc <- bcRecompute(bc, slot = "normalized")

## load nFeatures from seuv5
meta.data <- read.table("results/paad/seu_metadata_malignants.tsv")
bc@meta.data <- meta.data

bc <- bcRegressOut(bc = bc, vars.to.regress = c("nFeature_RNA"))

## export bcmat for seuv5
saveRDS(object = bc@normalized, file = "results/paad/bc_normalized_mat_rgrss.rds")

## sanity check
bc <- bcUMAP(bc = bc, pc = 10, k.neighbors = 20, npcs = 15, res = 1)
bcClusters(bc, idents = "bc_clusters_res.1")
