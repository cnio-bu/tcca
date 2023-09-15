library("Seurat")
library("tidyverse")

coding_genes <- data.table::fread(
    "/storage/scratch01/shared/projects/bc-meta/reference/hgnc_gene_with_protein_product_2023-03-22.tsv"
    ) %>%
    pull(symbol)


regenerate_chunk_with_coding <- function(x){
    
    x <- readRDS(x)
    ## Get rid of genes not in our annotation
    this_counts <- LayerData(x, layer = "counts")
    genes_to_keep <- intersect(rownames(this_counts), coding_genes)
    this_counts <- this_counts[genes_to_keep, ]
    
    # Filter columns to standarize
    meta.data <- x@meta.data %>%
        select(
            sample,
            percent.mt,
            percent.ribo,
            nFeature_RNA,
            nCount_RNA,
            Phase,
            G2M.Score,
            S.Score,
            cell_type
        )

    chunk_annotated <- Seurat::CreateSeuratObject(
        counts = this_counts,
        meta.data = meta.data,
        project = "pancancer"
    )
    
    return(chunk_annotated)
}

## Load all chunks
all_chunks <- list.files(
    "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/",
    full.names = TRUE,
    pattern = "all_objects_v4"
    )

all_chunks_annotated <- lapply(X = all_chunks, regenerate_chunk_with_coding)

full_object <- merge(
    x = all_chunks_annotated[[1]],
    y = all_chunks_annotated[2:length(all_chunks_annotated)]
)

saveRDS(object = full_object, file = "/storage/scratch01/shared/projects/bc-meta/single_cell/full_merged.rds")