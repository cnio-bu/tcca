#options(timeout = max(300, getOption("timeout")))

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("CellBench", "BiocStyle", "scater"))

# devtools::install_github("Oshlack/speckle")
# devtools::install_github("Oshlack/cellXY")


library(Seurat)
library(BPCells)
library(speckle)
library(cellXY)
library(CellBench)
library(org.Hs.eg.db)
library(dplyr)

seu <- readRDS("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl2/seu_lvl2.rds")
subset <- subset(seu, subset = sex == "")
joined <- JoinLayers(subset)

mat <- as.matrix(joined[["RNA"]]$counts)

sex <- classifySex(mat, genome = "Hs", qc = FALSE)

saveRDS(sex, "/storage/scratch01/users/mgonzalezb/bc-meta/sex_inference.Rds")

sex <- readRDS("/storage/scratch01/users/mgonzalezb/bc-meta/sex_inference.Rds")

sex <- sex %>%
bind_cols(subset@meta.data) %>%
mutate(study_patient = paste(study, patient, sep = "__")) %>%
mutate(study_patient = ifelse(patient == "" | patient == "ccl", paste(study, sample, sep = ), study_patient))

tab_pct <- t(apply(table(sex$study_patient, sex$prediction), 1, function(x) x/sum(x)))

sex_assign <- data.frame(
  study_patient = rownames(tab_pct),
  sex = apply(tab_pct, 1, function(row) {
    pred <- colnames(tab_pct)[which(row >= 0.8)]
    if (length(pred) > 0) pred else NA  # Assign NA if no sex meets the 80% criterion
  }),
  stringsAsFactors = FALSE
)

sex_full <- seu@meta.data %>% 
mutate(study_patient = ifelse(patient == "" | patient == "ccl", 
                                paste(study, sample, sep = "__"), 
                                paste(study, patient, sep = "__"))) %>%
left_join(sex_assign, by = "study_patient") %>%
mutate(sex = ifelse(sex.x == "", sex.y, sex.x)) %>%
mutate(sex = ifelse(sex == "Male" | sex == "m", "m", "f")) %>%
select(sex)

seu$sex <- as.vector(sex_full)
saveRDS(seu, "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl2/seu_lvl2_sex_infered.rds")