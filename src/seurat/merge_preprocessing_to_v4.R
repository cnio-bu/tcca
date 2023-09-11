library("Seurat")
library("tidyverse")

malignant_studies <- list.files(path = "raw/malignants/", full.names = TRUE)

## Load the very first study
first_study <- readRDS(malignant_studies[1])

## first merge
seu_first_sample <- first_study[[1]]

seu_first_study_merged <- merge(
  x = seu_first_sample,
  y = first_study[2:length(first_study)],
  merge.data = TRUE,
  )

## With the first study aggregated, start agg 1 study at a time.--
for(study in malignant_studies[2:length(malignant_studies)]){
  this_study <- readRDS(study)
  
  seu_first_study_merged <- merge(
    x = seu_first_study_merged,
    y = this_study
  )
}

## Save for later transform in v5
saveRDS(object = seu_first_study_merged, file = "all_objects_v4_merged.rds")

