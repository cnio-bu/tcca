library(BPCells)
library(tidyverse)
library(Seurat)

all_mats <- list.dirs(
    path = "/storage/scratch01/users/mgonzalezb/bc-meta/perception",
    full.names = TRUE
    )


# yay!
cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

all_mats <- all_mats[2:length(all_mats)]

mats <- map(
    all_mats,
    open_matrix_dir
)

full_mat <- do.call(cbind.fill, mats)
full_mat <- as(full_mat, "sparseMatrix")

write_matrix_dir(
    mat = full_mat,
    dir = "/storage/scratch01/users/mgonzalezb/bc-meta/perception/full_mat_beyondcell"
    )

rm(full_mat, mats)
gc()

mat <- open_matrix_dir(dir = "/storage/scratch01/users/mgonzalezb/bc-meta/perception/full_mat_beyondcell")
