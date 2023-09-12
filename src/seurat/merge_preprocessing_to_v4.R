library("Seurat")
library("tidyverse")

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/malignant")

malignant_studies <- list.files(path = "./", full.names = TRUE)

## Load the very first study
first_study <- readRDS(malignant_studies[1])

## first merge
print("Perform first merge")
seu_first_sample <- first_study[[1]]

seu_first_study_merged <- merge(
  x = seu_first_sample,
  y = first_study[2:length(first_study)],
  merge.data = TRUE,
  )

print("Load all objects")
all_studies <- do.call(readRDS, malignant_studies[2:length(malignant_studies)])
print("Flattening list")
all_studies <- unlist(all_studies, recursive = FALSE)


## With the first study aggregated, start agg 1 study at a time.--
#for(study in malignant_studies[2:length(malignant_studies)]){
#  this_study <- readRDS(study)
#  print(paste0("Merging ", study))  
#  seu_first_study_merged <- merge(
#    x = seu_first_study_merged,
#    y = this_study
#  )
#}


print("Performing merge...")
seu_first_study_merged <- merge(
  x = seu_first_study_merged,
  y = all_studies
)

## Save for later transform in v5
saveRDS(object = seu_first_study_merged, file = "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/all_objects_v4_merged.rds")
