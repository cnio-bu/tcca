library(Seurat)
library(BPCells)
setwd("/storage/scratch01/shared/projects/bc-meta/")

bc <- readRDS("beyondcell/beyondcell_pancancer.rds")

# Subset GBM samples from gbm_nourhan_abdelfattah study
bc_gbm <- subset(bc, subset = study == "gbm_nourhan_abdelfattah")

# Write h5ad object from BPCells bc matrix as input to SCellBow
write_matrix_anndata_hdf5(mat = bc_gbm[["RNA"]]$data, 
                          path = "single_cell/scellbow/beyondcell/gbm_bcs_target.h5ad")

# Save the metadata for those cells
write.table(bc_gbm@meta.data,
            file = "single_cell/scellbow/beyondcell/gbm_metadata_target.tsv",
            row.names = TRUE,
            sep = "\t")

# Subset expression data
seu <- readRDS("single_cell/seurat/v5/lvl2/seu_lvl2_sex_inferred.rds")
malignant <- subset(seu, subset = malignancy == TRUE)
colnames(malignant) <- paste0("c", c(1:ncol(malignant)))
write_matrix_anndata_hdf5(mat = malignant[["RNA"]]$counts.gbm_nourhan_abdelfattah,
                          path = "single_cell/scellbow/expression/gbm_expr_target.h5ad")

# Save the metadata for those cells
write.table(malignant@meta.data[malignant$study == "gbm_nourhan_abdelfattah", ],
            file = "single_cell/scellbow/expression/gbm_metadata_target.tsv",
            row.names = TRUE,
            sep = "\t")