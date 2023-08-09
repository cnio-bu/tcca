library("beyondcell")
library("Seurat")
library("tidyverse")

## SNAKEMAKE I/O
full_seurat_list <- snakemake@input[["seurat_list"]]
malignant_list <- snakemake@output[["malignant_list"]]
bc_list        <- snakemake@output[["bc_list"]]
report_file         <- snakemake@output[["report"]]

## Annotate the cell cycle score and cell phase for each cell, from a list
## of seurat objects
annotate_cell_cycle <- function(sc){
    sc <- CellCycleScoring(
        object = sc,
        s.features = cc.genes$s.genes,
        g2m.features = cc.genes$g2m.genes
    )
    return(sc)
}

## For this study, filter out non malignant cells
filter_malignant <- function(sc) {
  
  types_to_keep <- "Cancer Epithelial"
  if (sum(sc@meta.data$celltype_major == types_to_keep) > 1) {
    sc_filtered <- subset(x = sc, subset = celltype_major == types_to_keep)
    return(sc_filtered)
    
  } else {
    return(NULL)
  }
  
}

## Calculate bcscores sample wise, for each sample, for malignant pops. only
gs <- beyondcell::GetCollection(SSc, n.genes = 250, include.pathways = FALSE)

get_bcscores <- function(sc){
    bc <- bcScore(sc = sc, gs = gs, expr.thres = 0.1)
    
    # Do not allow NaNs
    bc@normalized[is.na(bc@normalized)] <- 0
    bc <- bcRecompute(bc, slot = "normalized")
    bc <- bcRegressOut(bc = bc, vars.to.regress = c("nFeature_RNA"))
    return(bc)
}

## Perform operations over list
seu <- readRDS(file = full_seurat_list)
seu <- lapply(X = seu, FUN = annotate_cell_cycle)

malignants <- lapply(X = seu, FUN = filter_malignant)

## Get rid of the NULL elements
malignants[sapply(malignants, is.null)] <- NULL
# Get rid of samples with < 100 malignant cells
malignants <- malignants[sapply(malignants, ncol) >= 100]

bcs <- lapply(X = malignants, FUN = get_bcscores)

## Save results objects
saveRDS(object = malignants, file = malignant_list)
saveRDS(object = bcs, file = bc_list)

## Generate and save reports
single_cell_report <- data.frame(
    sample = sapply(seu, FUN = function(x){unique(x@meta.data$orig.ident)}),
    cells = sapply(seu, FUN = function(x){ nrow(x@meta.data)})
)

bc_report <- data.frame(
    sample = sapply(malignants, FUN = function(x){ unique(x@meta.data$orig.ident)}),
    malignants = sapply(malignants, FUN = function(x){ nrow(x@meta.data)}),
    drug_sigs = sapply(bcs, FUN = function(x){ nrow(x@normalized)})
) 

single_cell_report <- single_cell_report %>%
    full_join(y = bc_report, by = "sample") %>%
    replace_na(list(malignants = 0, drug_sigs = 0))

write.table(
    x = single_cell_report,
    file = report_file,
    sep = "\t",
    row.names = FALSE
    )
