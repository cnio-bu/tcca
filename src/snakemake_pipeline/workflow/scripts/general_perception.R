library("Seurat")
library("Matrix")
library("BPCells")
library("caret")

## SNAKEMAKE I/O
malignant_list <- snakemake@input[["malignant_list"]]
drug_models           <- snakemake@input[["drug_models"]]
perception_mat <- snakemake@output[["perception_mat"]]


# Given a drug response model and an expression matrix, this function returns
# a viability matrix
viability_from_model <- function(infunc_DrugName,
                                 infunc_model,
                                 infunc_dataset) {
  coefnames <- infunc_model$coefnames
  overlap_genes <- intersect(coefnames, rownames(infunc_dataset))
  zero_genes <- setdiff(coefnames, rownames(infunc_dataset))
  zero_mat <- Matrix(
    0,
    nrow = length(zero_genes),
    ncol = ncol(infunc_dataset),
    sparse = TRUE
  )
  rownames(zero_mat) <- zero_genes
  
  # FOI stands for features of interest
  infunc_dataset_FOI <- infunc_dataset
  infunc_dataset_FOI <- rbind(infunc_dataset_FOI, zero_mat)
  infunc_dataset_FOI_t <- data.frame(t(infunc_dataset_FOI))
  
  # Predict viability score based on drug predictive model
  viability_score <- predict(infunc_model, infunc_dataset_FOI_t)
  return(viability_score)
}


## Perform operations over list
mals <- readRDS(file = malignant_list)
drug_models <- readRDS(drug_models)

# Compute predictive killing score per cell in each study separately
killing_mats <- lapply(1:length(mals), function(sample_n) {
  sample <- names(mals)[sample_n]
  mat <- mals[[sample_n]][["RNA"]]@data
  
  killing_mat_list <- lapply(1:length(drug_models), function(drug) {
    viability_from_model(
      infunc_DrugName = names(drug_models)[drug],
      infunc_model = drug_models[[drug]]$model,
      infunc_dataset = mat
    )
  })
  print(sample)
  killing_mat <- do.call(rbind, killing_mat_list)
  killing_mat <- as(killing_mat, "sparseMatrix")
  return(killing_mat)
})

# Transform matrices into Seurat v5 objects
names(killing_mats) <- names(mals)
full_killing_mat <- do.call(cbind, killing_mats)
full_killing_mat <- as(full_killing_mat, "sparseMatrix")

write_matrix_dir(mat = full_killing_mat,
                 dir = perception_mat,
                 overwrite = TRUE)
