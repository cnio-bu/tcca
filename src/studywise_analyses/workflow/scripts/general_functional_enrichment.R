library("beyondcell")
library("Seurat")

## SNAKEMAKE I/O
malignant_list <- snakemake@input[["malignant_cells"]]
gset           <- snakemake@input[["gsets"]]
bc_list        <- snakemake@output[["bc_list"]]

## Calculate bcscores sample wise, for each sample, for malignant pops. only
gs <- beyondcell::GenerateGenesets(x = gset, perform.reversal = FALSE)

## Safe way to calculate scores for datasets in which the data == count
## slot because of starting from norm.mat.
get_bcscores <- function(sc){
    mat <- as.matrix(GetAssayData(sc, slot = "data"))
    meta <- sc@meta.data
    bc <- bcScore(sc = mat, gs = gs, expr.thres = 0.1)
    bc@meta.data <- meta
    # Do not allow NaNs
    bc@normalized[is.na(bc@normalized)] <- 0
    bc <- bcRecompute(bc, slot = "normalized")
    bc <- bcRegressOut(bc = bc, vars.to.regress = c("nFeature_RNA"))
    return(bc)
}

## Perform operations over list
mals <- readRDS(file = malignant_list)
bcs <- lapply(X = mals, FUN = get_bcscores)

## Save results objects
saveRDS(object = bcs, file = bc_list)