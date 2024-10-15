library(beyondcell)

mat <- readRDS("results/paad/normalized_seu_malignant_mat.rds")

gs <- GetCollection(SSc, include.pathways = FALSE)
bc <- bcScore(sc = mat, gs = gs, expr.thres = 0.1)
bc@normalized[is.na(bc@normalized)] <- 0
bc <- bcRecompute(bc, slot = "normalized")

## load nFeatures from seuv5
meta.data <- read.table("results/paad/seu_metadata_malignants.tsv")
bc@meta.data <- meta.data

bc <- bcRegressOut(bc = bc, vars.to.regress = c("nFeature_RNA"))

## export bcmat for seuv5
saveRDS(object = bc@normalized, file = "results/paad/bc_normalized_mat_rgrss.rds")
