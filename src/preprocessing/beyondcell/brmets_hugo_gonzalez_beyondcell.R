library("beyondcell")
library("Seurat")

cfg <- read.delim(file = "config/basics.tsv") 
keys <- cfg$key
cfg <- as.character(cfg$value)
names(cfg) <- keys

data_folder <- cfg["working_dir"]
data_folder <- paste0(data_folder, "/single_cell")

## FUNC. def.
calculate_cell_cycle <- function(sc) {
    sc <- CellCycleScoring(
        sc,
        s.features = s.genes,
        g2m.features = g2m.genes,
        set.ident = FALSE,
        search = FALSE,
        seed = 120394,
        verbose = FALSE
    )
    return(sc)
}

calculate_bc_scores <- function(sc) {
    bc <- bcScore(sc, gs = gs, expr.thres = 0.1)
    bc <- bcRegressOut(bc = bc, vars.to.regress = c("nFeature_RNA", "G2M.Score"))
    return(bc)
}


## Get Seurat gene sets for cell cycle.
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

## Load Seurat objs and genesets
malignant_sc <- readRDS(
    paste0(data_folder, "/obj/brmets_hugo_gonzalez/all_malignant.rds")
)

gs <- beyondcell::GenerateGenesets(
    paste0(cfg["working_dir"], "/reference/bc_signatures_cutoff.gmt"),
    include.pathways = FALSE,
    n.genes = 250,
    mode = c("up", "down")
)

## Remove samples with <= 50 cells. 
malignant_sc <- malignant_sc[sapply(malignant_sc, ncol) > 50]

## Calculate cell cycle scoring for all samples
malignant_sc <- lapply(malignant_sc, calculate_cell_cycle)

 # Calculate bc scores
bcs <- lapply(malignant_sc, calculate_bc_scores)

rm(malignant_sc)
gc()

## Prepare bc report
report <- data.frame(
    sample = names(bcs),
    cells = sapply(bcs, FUN = function(x){ ncol(x@normalized)}),
    sigs = sapply(bcs, FUN = function(x){ nrow(x@normalized)})
)


## If any sample has <= 450 sigs. we'll have to remove it
samples_to_keep <- report[report$sigs >= 450, "sample"]

bcs_filtered <- bcs[samples_to_keep]

## save everything
save_to <- paste0(cfg["working_dir"], "/beyondcell")
saveRDS(
    object = bcs_filtered,
    file = paste0(save_to, "/obj/brmets_hugo_gonzalez/bcs.rds")
)

write.table(
    x = report[samples_to_keep, ],
    file = paste0(save_to, "/obj/brmets_hugo_gonzalez/report.tsv"),
    sep = "\t",
    row.names = FALSE
)
