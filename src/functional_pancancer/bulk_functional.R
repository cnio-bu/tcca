library(edgeR)
library(GSVA)
library(tidyverse)

mat <- read.table("results/functional/pancancer_pseudobulk.tsv")
mat <- as.matrix(mat)

dge <- DGEList(counts = mat)
dge <- calcNormFactors(object = dge, method = "TMM")

## get norm. mat
mat_cpm <- cpm(
    dge,
    normalized.lib.sizes = TRUE,
    log = TRUE,
    prior.count = 1
)

dim(mat_cpm)

## gsva
gsets <- GSEABase::getGmt(con = "reference/combined_gsets_functional.gmt")

gsets_fix <- list()
for(gset in gsets){
    gset@geneIds <- gset@geneIds[gset@geneIds != ""]
    gsets_fix <- c(gsets_fix, gset)
}

gsets_fixed <- GSEABase::GeneSetCollection(gsets_fix)

## gsva parameters
gsvapar <- gsvaParam(
    exprData = mat_cpm,
    geneSets = gsets_fixed,
    kcdf = "Gaussian",
    maxDiff = TRUE
    )

gsva_es <- gsva(gsvapar)

## test
ComplexHeatmap::Heatmap(
    matrix = gsva_es,
    show_column_names = FALSE,
    clustering_distance_rows = "pearson",
    clustering_distance_columns = "pearson",
    row_split = 10,
    column_split = 10,
    row_names_gp = gpar(fontsize = 7)
    )


