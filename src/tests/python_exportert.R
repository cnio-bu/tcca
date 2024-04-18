library(BPCells)
library(tidyverse)

mat <- BPCells::open_matrix_dir(dir = "results/lvl1/pancancer_merged_mat")
mat_to_ints <- BPCells::convert_matrix_type(mat, type = "uint32_t")

BPCells::write_matrix_anndata_hdf5(
    mat = mat_to_ints,
    path = "/raid/lsagarcia/seu_lvl1_counts.h5ad",
    # group = "counts"
)


BPCells::write_matrix_hdf5(
    mat = mat_to_ints,
    path = "/raid/lsagarcia/seu_lvl1_counts.hdf5",
    group = "counts",
    compress = TRUE
)

## export the metadata as is
seu <- readRDS("results/lvl1/seu_lvl1_merged.Rds")

meta.data <- seu@meta.data %>%
    rownames_to_column("rowname_cell")

write_tsv(x = meta.data, file = "/raid/lsagarcia/metadata_lvl1_h5.tsv")
