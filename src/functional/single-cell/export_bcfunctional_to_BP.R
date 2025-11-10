library(BPCells)
library(beyondcell)
library(Matrix)
library(tidyverse)


# yay!
cbind.fill <- function(...){
    nm <- list(...) 
    nm<-lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

all.studies <- list.files(
    path = "/storage/scratch01/shared/projects/bc-meta/single_cell/functional",
    pattern = "*.rds",
    full.names = TRUE
    )

all.studies <- all.studies[c(1:5, 7:37)]

for(study in all.studies){

    study_name <- basename(study)
    study_name <- gsub(x=study_name, pattern=".rds", replacement="")
    samples <- readRDS(study)
    print(study_name)
    study_matrices <- list()
    study_meta <- list()

    for(sample in samples){
        sample_name <- unique(sample@meta.data$sample)
        compressed_mat <- as(sample@normalized, "sparseMatrix")
        study_matrices <- append(study_matrices, compressed_mat)
        study_meta <- append(study_meta, list(sample@meta.data))
    }
    if(length(study_matrices) > 1){
        full_mat <- do.call(cbind.fill, study_matrices)
        full_mat <- as(full_mat, "sparseMatrix")
        full_meta <- bind_rows(study_meta)
     }else{
        full_mat <- study_matrices[[1]]
        full_mat <- as(full_mat, "sparseMatrix")
        full_meta <- study_meta[[1]]
}

    write_matrix_dir(
        mat = full_mat,
        dir = paste0("/storage/scratch01/shared/projects/bc-meta/functional/", study_name),
        overwrite=TRUE
        )
    
    write.table(x = full_meta, file = paste0(
        "/storage/scratch01/shared/projects/bc-meta/functional/",
        study_name,
        ".tsv"
        ),
    )
}
