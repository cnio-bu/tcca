#!/usr/bin/env Rscript

#sbatch -c 20 --job-name=scthx -o log_scthx.txt -e error_scthx.txt --mem=80G -t 1000 --wrap "Rscript /home/lserranor/bc-meta/src/sctherapy/sctherapy_fullcohort.R"

library(BPCells)
library(Seurat)
library(tidyverse)
library(httr)
library(logger)
library(jsonlite)

invisible(lapply(c("https://raw.githubusercontent.com/kris-nader/scTherapy/main/R/identify_healthy_mal_v5.R",
                   "https://raw.githubusercontent.com/kris-nader/scTherapy/main/R/identify_subclones_v5.R",
                   "https://raw.githubusercontent.com/kris-nader/scTherapy/main/R/predict_compounds.R"),source))


setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl2/")
where_to_save <- "/storage/scratch01/shared/projects/bc-meta/single_cell/sctherapy/"

if (!dir.exists(where_to_save)) {
  dir.create(where_to_save)
}

##### FUNCTIONS

extract_seu <- function(study) {
  path <- paste0("./", study, "_v5/")
  mat <- open_matrix_dir(path)
  mat <- mat[, colnames(mat) %in% rownames(full_metadata)]
  seu <- CreateSeuratObject(mat,
                            meta.data = full_metadata)
  return(seu)
}

make_sample_list <- function (seu) {
  seu_list <- list()
  for (sam in unique(seu$sample)){
    sub_seu <- subset(seu, sample == sam)
    seu_list <- c(seu_list, list(sub_seu))
  }
  names(seu_list) <- unique(seu$sample)
  return(seu_list)
}

gene_list <- "https://raw.githubusercontent.com/kris-nader/scTherapy/main/geneinfo_beta_input.txt"
gene_info <- data.table::fread(gene_list) %>% as.data.frame()

get_malignant_deg <- function(patient_sample){
  
  # Normalize data
  patient_sample <- NormalizeData(patient_sample, normalization.method = "LogNormalize", scale.factor = 10000)
  patient_sample <- FindVariableFeatures(patient_sample, selection.method = "vst", nfeatures = 2000)
  
  ## Get DEG
  patient_sample@meta.data <- patient_sample@meta.data %>%
    mutate(ensemble_output = ifelse(malignancy, "malignant", "healthy"))
  
  Idents(patient_sample) <- patient_sample@meta.data$ensemble_output
  
  deg_malignant <- FindMarkers(object = patient_sample, 
                               ident.1 = "malignant", 
                               ident.2 = "healthy", 
                               min.pct = -Inf, logfc.threshold = -Inf,
                               min.cells.feature = 1, min.cells.group = 1,
                               test.use="wilcox")
  
  ## Filter DEG					 
  deg_malignant <- deg_malignant %>%
    mutate(gene_symbol = rownames(.)) %>% inner_join(gene_info, by = "gene_symbol") %>%
    filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
  
  deg_malignant <- setNames(as.list(deg_malignant$avg_log2FC), deg_malignant$gene_symbol)      
  
  return(deg_malignant)
}

get_subclones_deg <- function(patient_sample){
  
  # Normalize data
  patient_sample <- NormalizeData(patient_sample, normalization.method = "LogNormalize", scale.factor = 10000)
  patient_sample <- FindVariableFeatures(patient_sample, selection.method = "vst", nfeatures = 2000)
  
  # Get DEG for all subclones in the sample
  all_subclones <- setdiff(unique(patient_sample$ensemble_output), "healthy")
  
  deg_all_subclones <- list()
  
  for (subclone in all_subclones){
    
    ## Get DEG
    Idents(patient_sample) <- patient_sample@meta.data$ensemble_output
    
    deg_subclone <- FindMarkers(object = patient_sample, 
                                ident.1 = subclone, 
                                ident.2 = "healthy", 
                                min.pct = -Inf, logfc.threshold = -Inf,
                                min.cells.feature = 1, min.cells.group = 1,
                                test.use="wilcox")
    
    ## Filter DEG					 
    deg_subclone <- deg_subclone %>%
      mutate(gene_symbol = rownames(.)) %>% inner_join(gene_info, by = "gene_symbol") %>%
      filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
    
    deg_subclone <- setNames(as.list(deg_subclone$avg_log2FC), deg_subclone$gene_symbol)
    
    deg_all_subclones <- c(deg_all_subclones, list(deg_subclone))
  }
  
  names(deg_all_subclones) <- all_subclones
  
  return(deg_all_subclones)
}


##### FULL-COHORT METADATA
## Generate metadata + clonality info for the full cohort
full_seu <- readRDS("seu_lvl2_sex_inferred.rds")
metadata <- full_seu@meta.data %>%
  rownames_to_column("original_barcode")

clonality <- read.table("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata/full_clonality_table_lvl2.tsv", row.names = NULL)
duplicated <- clonality[clonality$original_barcode %in% clonality$original_barcode[duplicated(clonality$original_barcode)], ]
clonality <- clonality %>%
  filter(!(original_barcode %in% duplicated$original_barcode)) %>% # Remove 16 cells with duplicated names
  dplyr::select(original_barcode, scevan_prediction, confidentNormal, scevan_subclone) 

full_metadata <- left_join(metadata, clonality, by = "original_barcode")

## Keep only cells in agreement for better quality of the prediction
full_metadata <- full_metadata %>%
  filter(malignancy == F & scevan_prediction == "normal" | malignancy == T & scevan_prediction == "tumor") 

full_metadata <- full_metadata %>%
  mutate(subclone_name = paste(study, sample, scevan_subclone, sep = "."),
         ensemble_output = ifelse(malignancy == FALSE, "healthy", subclone_name)) %>%
  column_to_rownames(var = "original_barcode")


##### CODE
all_study_tables <- list()

for (study_name in unique(full_metadata$study)) {
  ## Get Seurat object
  study_seu <- extract_seu(study_name) 
  
  ## Make Seurat per sample list
  study_seu <- make_sample_list(study_seu)
  
  ## Remove samples with no malignant cells
  study_seu <- study_seu[sapply(study_seu, function(x) sum(x$malignancy == FALSE) > 0)]
  
  ## If all samples were removed, move on to next study
  if (length(study_seu) == 0) {
    next
  }

  ## Get DEG lists and save them
  study_deg <- lapply(study_seu, get_subclones_deg)
  
  saveRDS(study_deg, paste0(where_to_save, study_name, ".rds"))
  
  ## Predict drug:dose combinations
  error_log_file <- paste0(where_to_save, "drug_prediction_error_log.txt")

  study_drugs <- lapply(study_deg, function(sample) {
    setNames(lapply(names(sample), function(subclone_name) {
      tryCatch({
        predict_drugs(degs_list = sample[[subclone_name]], exclude_low_confidence = TRUE)
      }, error = function(e) {
        message <- paste(Sys.time(), "Error en subclone", subclone_name, ":", e$message, "\n")
        cat(message, file = error_log_file, append = TRUE)
        return(NULL)  ## Return NULL in case of error
      })
    }), names(sample))  
  })

  ## Generate un-nested table, save individually and add to common list

  final_dataframe <- do.call(rbind, lapply(names(study_drugs), function(sample_name) {
    valid_subclones <- Filter(Negate(is.null), study_drugs[[sample_name]])  ## Remove NULLs
    do.call(rbind, lapply(names(valid_subclones), function(subclone_name) {
      df <- valid_subclones[[subclone_name]]
      df$Sample <- sample_name
      df$Subclone <- subclone_name
      return(df)
    }))
  }))
  
  write.table(final_dataframe, paste0(where_to_save, study_name, ".tsv"), sep = "\t")
  
  all_study_tables <- c(all_study_tables, list(final_dataframe))
}

full_table <- do.call(rbind, all_study_tables)

write.table(full_table, paste0(where_to_save, "full_table_drug_prediction.tsv"), sep = "\t")
