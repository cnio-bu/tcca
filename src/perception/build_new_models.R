working_directory="/storage/scratch01/users/mgonzalezb/bc-meta/perception_newmodeles/"
setwd(working_directory)

onTarget <- readRDS('Data/DepMapv12.RDS')
# Function to load a RMD file to another
ksource <- function(x, ...) {
  library(knitr)
  source(purl(x, output = tempfile()), ...)
}
ksource('./Tools/step0A_functions_needed.Rmd')
ksource('./Tools/step0B_functions_needed.Rmd')
genesUsed_toBuild=readRDS('Data/genesUsed_toBuild.RDS')
genes_across_scRNA_datasets_ofInterest=readRDS('Data/genes_across_scRNA_datasets_ofInterest.RDS')

seu <- readRDS("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5/lvl1/seu_lvl1.rds")
DOI <- onTarget$secondary_screen_drugAnnotation$CommonName
cancerType='PanCan'
# What features are we using?
genesUsed_toBuild=Reduce(intersect,list(rownames(onTarget$expression_20Q4),
                                        rownames(onTarget$scRNA_complete),
                                        rownames(seu)))

# How many features should the model use?
Possible_Number_of_features=nrow(onTarget$expression_rnorm) * seq(0.0005, 0.02, length = 5)

features_bulk_all_cancer_type_Cavg <- readRDS("features_bulk_all_cancer_type_Cavg.Rds")
Raw_models_output <- list();
Performance_by_features <- list(); for_output_lung_Test_vglm <- list()
# Remove drugs with NA
names(features_bulk_all_cancer_type_Cavg) <- DOI
features_bulk_all_cancer_type_Cavg[is.na(features_bulk_all_cancer_type_Cavg)] <- NULL
DOI <- names(features_bulk_all_cancer_type_Cavg)
cl <- parallel::makeForkCluster(10)
doParallel::registerDoParallel(cl)

# Iterate over the DOI (drugs of interest)
for (i in 1:length(DOI)) {
  # Initialize variables
  Raw_models_output <- list()
  id_counter <- 1
  # For each given drug; iterate over a parameter of 'features count'
  for (infunc_k_features_grid in Possible_Number_of_features ) {
    # Build the model
    Raw_models_output[[id_counter]] <- build_on_BULK_v2(
      infunc_drugName = DOI[i],
      infunc_cancerType = cancerType,
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
  print(sapply(Raw_models_output, function(x) x$performance_in_scRNA))
  Performance_by_features[[i]] <- lapply(Raw_models_output, function(x) x$performance_in_scRNA)

  # Store the Tuned model
  for_output_lung_Test_vglm[[i]] <- err_handle(Raw_models_output[[which.max(sapply(Raw_models_output, function(x) x$performance_in_scRNA)[2,])]])
  
  # Remove temporary variables and garbage collect
  rm(Raw_models_output)
  gc()
print(i)
}

# Stop parallel processing
parallel::stopCluster(cl)

# Assign names to the TUned models and save them
names(for_output_lung_Test_vglm) <- DOI
Tuned_models_output <- for_output_lung_Test_vglm
saveRDS(Tuned_models_output, paste('Data/new_models_DOI', Sys.time(), '.rds',sep=''))