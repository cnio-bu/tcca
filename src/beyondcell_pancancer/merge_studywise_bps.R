library(BPCells)
library(tidyverse)

all_mats <- list.dirs(
    path = "results/beyondcell_bp",
    full.names = TRUE
    )


# yay!
cbind.fill<-function(...){
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}


all_mats <- all_mats[2:37]

mats <- map(
    all_mats,
    open_matrix_dir
)

full_mat <- do.call(cbind.fill, mats)
full_mat[is.na(full_mat)] <- 0
full_mat <- as(full_mat, "sparseMatrix")

pcs <- prcomp(x = full_mat, center = TRUE, scale. = TRUE)
