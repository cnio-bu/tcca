library(Seurat)
library(BPCells)


setwd("/home/mgonzalezb/bc-meta/PERCEPTION/")


## STEP1: loading neccesary files
# Function to load a RMD file to another
onTarget <- readRDS('Data/DepMapv12.RDS')
ksource <- function(x, ...) {
  library(knitr)
  source(purl(x, output = tempfile()), ...)
}
ksource('Tools/step0A_functions_needed.Rmd')
ksource('Tools/step0B_functions_needed.Rmd')
genesUsed_toBuild = readRDS('Data/genesUsed_toBuild.RDS')
genes_across_scRNA_datasets_ofInterest = readRDS('Data/genes_across_scRNA_datasets_ofInterest.RDS')

# Load already built models for all FDA approved drugs
approved_drugs_model <- readRDS('Data/FDA_approved_drugs_models.RDS')

# Load level 1 seurat object with all malignant cells from TCCA
seu.lvl1 <- readRDS("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl1/seu_lvl1.rds")

# Given a drug response model and an expression matrix, this function returns a viability matrix
viability_from_model <- function(infunc_DrugName,
                                 infunc_model,
                                 infunc_dataset){
    
    coefnames <- infunc_model$coefnames
    overlap_genes <- intersect(coefnames, rownames(infunc_dataset))
    # zero_genes <- setdiff(coefnames, rownames(infunc_dataset))
    # zero_mat <- matrix(0, nrow = length(zero_genes), ncol = ncol(infunc_dataset))
    # rownames(zero_mat) <- zero_genes

    # FOI stands for features of interest
    infunc_dataset_FOI <- infunc_dataset
    # infunc_dataset_FOI <- rbind(infunc_dataset_FOI, zero_mat)
    infunc_dataset_FOI_t <- data.frame(t(infunc_dataset_FOI))

    # Predict viability score based on drug predictive model
    viability_score <- predict(infunc_model,
                          infunc_dataset_FOI_t)
    return(viability_score)
}


## STEP2: Build PERCEPTION models
# Drug of interest to build the models (include some drugs from SSc collection in Beyondcell)
load("SSc.RData")
ssc_bc <- tolower(unique(SSc@info$preferred.drug.names))
doi <- unique(c(intersect(onTarget$drugCategory$name, ssc_bc), names(approved_drugs_model)))

# What features are we using?
genesUsed_toBuild <- Reduce(intersect,list(rownames(onTarget$expression_20Q4),
                                        rownames(onTarget$scRNA_complete),
                                        rownames(seu.lvl1)))

# Normalization of bulk and sc-expression
onTarget$scRNA_subset_rnorm <- rank_normalization_mat(onTarget$scRNA_complete[genesUsed_toBuild,])
onTarget$expression_rnorm <- rank_normalization_mat(onTarget$expression_20Q4[genesUsed_toBuild,])

# How many features should the model use?
possible_number_of_features <- nrow(onTarget$expression_rnorm) * seq(0.0005, 0.02, length = 5)

# Feature Selection/Ranking using Bulk
features_bulk_all_cancer_type_Cavg <- run_parallel_feature_ranking_bulk(id_cancerType = 'PanCan',
                                                                        infunc_DrugsToUse = doi,
                                                                        infunc_GOI=genesUsed_toBuild,
                                                                        resp_measure_mode = 'AUC',
                                                                        infunc_exclude_cancer = 'PanCan')

# How feature ranking looks like?
head(features_bulk_all_cancer_type_Cavg[[1]], 5)
names(features_bulk_all_cancer_type_Cavg) <- doi
na_drug_models <- sapply(features_bulk_all_cancer_type_Cavg, function(x) (length(x) == 1 && is.na(x)))
features_bulk_all_cancer_type_Cavg <- features_bulk_all_cancer_type_Cavg[!na_drug_models]

# Initialize list to save perception output. Provide cores for parallelization 
doi <- names(features_bulk_all_cancer_type_Cavg)
raw_models_output <- list()
performance_by_features <- list()
for_output_tcca <- list()
cl <- parallel::makeForkCluster(4)
doParallel::registerDoParallel(cl)


# Build and tune PERCEPTION
# Iterate over the doi (drugs of interest)
for (i in 1:length(doi)) {
  # Initialize variables
  raw_models_output <- list()
  id_counter <- 1
  # For each given drug; iterate over a parameter of 'features count'
  for (infunc_k_features_grid in possible_number_of_features) {
    # Build the model
    raw_models_output[[id_counter]] <- build_on_BULK_v2(
      infunc_drugName = doi[i],
      infunc_cancerType = 'PanCancer',
      infunc_features = rownames(features_bulk_all_cancer_type_Cavg[[i]]),
      single_best = rownames(features_bulk_all_cancer_type_Cavg[[i]])[1],
      k_features = infunc_k_features_grid,
      mode = 'AUC',
      model_type = 'glmnet',
      exclude_cancer = 'PanCan'
    )
    id_counter <- id_counter + 1
  }

  # Print and store the performance
  print(sapply(raw_models_output, function(x) x$performance_in_scRNA))
  performance_by_features[[i]] <- lapply(raw_models_output, function(x) x$performance_in_scRNA)

  # Store the Tuned model
  for_output_tcca[[i]] <- err_handle(raw_models_output[[which.max(sapply(raw_models_output, function(x) x$performance_in_scRNA)[2,])]])
  
  # Remove temporary variables and garbage collect
  rm(raw_models_output)
  gc()
}
# Stop parallel processing
parallel::stopCluster(cl)

# Assign names to the Tuned models and save them
names(for_output_tcca) <- doi
na_drug_models <- sapply(for_output_tcca, function(x) (length(x) == 1 && is.na(x)))
for_output_tcca <- for_output_tcca[!na_drug_models]
tuned_models_output <- for_output_tcca
saveRDS(tuned_models_output, paste('./Data/PERCEPTION_models_doi_tcca', Sys.time(), '.RDS', sep = ''))

tuned_models_output <- readRDS('Data/PERCEPTION_models_doi_tcca2024-11-05 22:46:39.744596.RDS')
doi <- names(tuned_models_output)

# Step3: Model Performance to select drug predictive models
# Obtain the drug response models
drug_models <- tuned_models_output
corr_matrix <- cbind(
  bulk = unlist(lapply(drug_models,function(x)x$performance_in_bulk[2])),
  pseudo_bulk = unlist(lapply(drug_models,function(x)x$performance_in_pseudo_bulk[2])),
  sc = unlist(lapply(drug_models,function(x)x$performance_in_scRNA[2])))
rownames(corr_matrix) <- names(drug_models)

pval_matrix <- cbind(
  bulk = unlist(lapply(drug_models,function(x)x$performance_in_bulk[1])),
  pseudo_bulk = unlist(lapply(drug_models,function(x)x$performance_in_pseudo_bulk[1])),
  sc = unlist(lapply(drug_models,function(x)x$performance_in_scRNA[1]))
)
rownames(pval_matrix) <-  names(drug_models)


# Number of cell lines with AUC values used in validation
totalcellLines_used_in_validation <- unlist(lapply(
  drug_models, function(x)length(na.omit(x$predVSgroundTruth$pred_gt_mscRNA$Observed))))

# Select the significant drug models -->
significant_drugs_names <- names(which(corr_matrix[,3]>0.3 & pval_matrix[,3]<0.05)) 

# Drug models from PERCEPTION
significant_drugs_models <- drug_models[significant_drugs_names]


seu.lvl1 <- JoinLayers(seu.lvl1)

# STEP4: Compute predictive killing score per cell in each study separately
for (study_name in unique(seu.lvl1$study)){
    print(study_name)
    seu <- subset(seu.lvl1, subset = study == study_name)
    seu <- NormalizeData(seu,
                        normalization.method = "LogNormalize",
                        scale.factor = 10000)
    mat <- as.matrix(seu[["RNA"]]$data)
    # Predict viability from pre-built models
    killing_mat_list <- lapply(1:length(significant_drugs_names), function(x) {
    viability_from_model(
        infunc_DrugName = significant_drugs_names[x],
        infunc_model = significant_drugs_models[[x]]$model,
        infunc_dataset = mat)

})
    killing_mat <- do.call(rbind, killing_mat_list)
    killing_mat <- as(killing_mat, "sparseMatrix")
    write_matrix_dir(
        mat = killing_mat,
        dir = paste0("/storage/scratch01/users/mgonzalezb/bc-meta/perception/", study_name),
        overwrite = TRUE
        )
    
}