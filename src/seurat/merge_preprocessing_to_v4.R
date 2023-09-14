library("Seurat")
library("tidyverse")

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/malignant")

malignant_studies <- list.files(path = "./", full.names = TRUE)

print("Load all objects")
all_studies <- malignant_studies[1:length(malignant_studies)] %>%
  map(readRDS)

print("Flattening list")
all_samples <- unlist(all_studies, recursive = FALSE)

## Split in 100 samples chunks
all_chunks <- 1:length(all_samples)

print("Split chunks")
chunk_list <- split(all_chunks, ceiling(seq_along(all_chunks) / 100))
for(chunk in chunk_list){
  first_seu <- all_samples[[chunk[[1]]]]
  print(paste0("Merge from", " ", unique(first_seu$sample)))

  merged_seu <- merge(
    x = first_seu,
    y = all_samples[chunk[[2]]:chunk[length(chunk)]])

    print("Saving chunk")
    saveRDS(
      object =merged_seu,
      file = paste0(
        "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/all_objects_v4_merged_",
        max(chunk),
        ".rds"
        )
    )
}

print("DONE!")