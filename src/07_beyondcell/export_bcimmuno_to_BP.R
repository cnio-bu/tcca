library(BPCells)
library(Matrix)
library(tidyverse)

all.studies <- list.files(
    path = "/storage/scratch01/shared/projects/bc-meta/single_cell/beyondcell_immuno",
    pattern = "*.rds",
    full.names=TRUE
    )


for(study in all.studies){

    study_name <- basename(study)
    study_name <- gsub(x=study_name, pattern=".rds", replacement="")
    samples <- readRDS(study)
    print(study_name)
    study_matrices <- list()
    study_meta <- list()

    for (sample in samples) {
        sample_name <- unique(sample@meta.data$sample)
        compressed_mat <- as(sample@normalized, "sparseMatrix")
        study_matrices <- append(study_matrices, compressed_mat)
        study_meta <- append(study_meta, list(sample@meta.data))
    }
    if (length(study_matrices) > 1) {
        all_drugs <- unique(unlist(lapply(study_matrices, rownames)))
        study_matrices_aligned <- lapply(study_matrices, function(mat) {
            mat2 <- Matrix::Matrix(0,
                nrow = length(all_drugs),
                ncol = ncol(mat), 
                sparse = TRUE
            )
            rownames(mat2) <- all_drugs
            colnames(mat2) <- colnames(mat)
            mat2[rownames(mat), ] <- mat
            return(mat2)
        })
        full_mat <- do.call(cbind, study_matrices_aligned)
        full_mat <- as(full_mat, "sparseMatrix")
        full_meta <- bind_rows(study_meta)
        print(dim(full_mat))
    } else {
        full_mat <- study_matrices[[1]]
        full_mat <- as(full_mat, "sparseMatrix")
        full_meta <- study_meta[[1]]
    }

    write_matrix_dir(
        mat = full_mat,
        dir = paste0("/storage/scratch01/shared/projects/bc-meta/beyondcell_immuno/studywise_bpcells/", study_name),
        overwrite=TRUE
        )
    
    write.table(x = full_meta, file = paste0(
        "/storage/scratch01/shared/projects/bc-meta/beyondcell_immuno/studywise_bpcells/",
        study_name,
        ".tsv"
        ),
    )
    
}
