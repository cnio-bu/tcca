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