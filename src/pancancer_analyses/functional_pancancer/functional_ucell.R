library(Seurat)
library(BPCells)
library(UCell)
library(qusage)

# Read the gmt file with the 6 drug resistance mechanisms gene sets.
gmt_list <- read.gmt("/storage/scratch01/shared/projects/bc-meta/reference/combined_gsets_functional.gmt")

# Join bidirectional signatures in the same vector adding '+' to the genes in the upregulared gene set 
# and '-' to the genes in the downregulated gene set.
sig_names <- sub("_(UP|DOWN|DN)$", "", names(gmt_list))
bisigs <- sig_names[duplicated(sig_names)]

for (sig in unique(sig_names)){
    up <- paste0(sig, "_UP")
    down <-  if (paste0(sig, "_DOWN") %in% names(gmt_list)) paste0(sig, "_DOWN") else paste0(sig, "_DN")

    if (sig %in% bisigs){
        gmt_list[[up]] <- paste0(gmt_list[[up]], "+")
        gmt_list[[down]] <- paste0(gmt_list[[down]], "-")
        gmt_list[[sig]] <- c(gmt_list[[up]], gmt_list[[down]])
        gmt_list[[up]] <- NULL
        gmt_list[[down]] <- NULL

    }else{
        if(up %in% names(gmt_list)){
            gmt_list[[up]] <- paste0(gmt_list[[up]], "+")
            gmt_list[[sig]] <- gmt_list[[up]]
            gmt_list[[up]] <- NULL
        }

        if(down %in% names(gmt_list)){
            gmt_list[[down]] <- paste0(gmt_list[[down]], "-")
            gmt_list[[sig]] <- gmt_list[[down]]
            gmt_list[[down]] <- NULL
        }
    }
}

# Load level 2 Seurat object
seu <- readRDS("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl2/seu_lvl2.rds")
malignant <- subset(seu, subset = malignancy == TRUE)
malignant <- JoinLayers(malignant)
colnames(malignant) <- paste0("c", c(1:ncol(malignant)))

# Load the beyondcell pancancer object to select only cells where BCS was computed
bc <- readRDS("/storage/scratch01/users/mgonzalezb/bc-meta/beyondcell/results/beyondcell_pancancer_final_res.Rds")
cells <- colnames(bc)

# Subset the malignant cells with BCS score from the counts slot
malignant <- subset(malignant, cells = cells)

# Split malignant matrix into 3 matrices for computing UCell scores
full_mat <- malignant[["RNA"]]$counts
nrow(full_mat)/3
seu_ucell <- AddModuleScore_UCell(malignant, features = gmt_list)

scores_ucell <- seu_ucell@meta.data[, names(gmt_list)]

write.table(scores_ucell, "/storage/scratch01/users/mgonzalezb/bc-meta/functional/functional_ucell.tsv", sep = "\t")