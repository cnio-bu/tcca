library("BPCells")
library("Seurat")
library("tidyverse")

options(future.globals.maxSize = 1e9)
setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/raw")

coding_genes <- data.table::fread(
    "/storage/scratch01/shared/projects/bc-meta/reference/hgnc_gene_with_protein_product_2023-03-22.tsv"
    ) %>%
    pull(symbol)

## Get all files
malignant_studies <- list.files(path = "./", full.names = TRUE)

for(study_name in malignant_studies){
    options(Seurat.object.assay.version = "v4")

    seuv4_list <- readRDS(study_name)
    this_study <- basename(study_name)
    this_study <- gsub(pattern = ".rds", replacement = "", x = this_study)
    print(this_study)
    if (length(seuv4_list) > 1){
    	seu_merged <- merge(seuv4_list[[1]], y = seuv4_list[2:length(seuv4_list)])
	}else{seu_merged <- seuv4_list[[1]]
}
    rm(seuv4_list)
    gc()
    print("Merged")
    meta.data <- seu_merged@meta.data %>%
    select(any_of(c(
            "sample",
            "percent.mt",
            "percent.ribo",
            "nFeature_RNA",
            "nCount_RNA",
            "cell_type",
            "malignancy"
    ))) %>%
    as.data.frame()

    counts <- seu_merged[["RNA"]]$counts
    genes_to_keep <- intersect(rownames(counts), coding_genes)

    counts <- counts[genes_to_keep, ]

    rm(seu_merged)
    gc()

    options(Seurat.object.assay.version = "v5")

    write_matrix_dir(
        mat = counts,
        dir = paste0(
            "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/all_cell_types/",
            this_study,
            "_v5" 
        )
    )

    rm(counts)
    gc()

    counts.mat <- open_matrix_dir(paste0(
            "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/all_cell_types/",
            this_study,
            "_v5" 
        ))

    seu_v5 <- Seurat::CreateSeuratObject(
        counts = counts.mat,
        meta.data = meta.data,
        project = this_study
        )
    print("Seurat v5")
    saveRDS(
        object = seu_v5,
        file = paste0(this_study, "_v5.Rds"),
        destdir = "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/all_cell_types/" 
    )
}
