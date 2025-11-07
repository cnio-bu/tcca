library(Seurat)
library(BPCells)
library(tidyverse)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/")

summarize_cells <- function(metadata, level) {
  metadata %>%
    group_by(sample) %>%
    summarise(
      !!paste0("ncells.", level) := n(),
      !!paste0("ncells.mal.", level) := sum(malignancy == TRUE, na.rm = TRUE),
      !!paste0("ncells.tme.", level) := sum(malignancy == FALSE, na.rm = TRUE)
    )
}

studies <- list.dirs("./all_cell_types/")
studies <- grep("_v5", studies, value = TRUE)

# Remove study 'bmets_youmna_kfoury' that do not have any sample with more than 
# 100 malignant cells.
studies <- studies[!grepl("bmets_youmna_kfoury", studies)]

clinical <- read.table(
  "./v5/clinical_metadata_v4_clean.tsv",
  sep = "\t",
  header = TRUE
)
clinical <- as.data.frame(clinical) %>%
  filter(study != "brca_bhupinder_pal" |
           tumor_subtype != "predicted_tumour")

data.list_lvl1 <- c()
metadata.list_lvl1 <- c()
data.list_lvl2 <- c()
metadata.list_lvl2 <- c()
data.list_lvl3 <- c()
metadata.list_lvl3 <- c()

report_df <- data.frame(
  study = character(),
  sample = character(),
  ncells.mal.lvl3 = numeric(),
  ncells.tme.lvl3 = numeric(),
  ncells.lvl3 = numeric(),
  ncells.mal.lvl2 = numeric(),
  ncells.tme.lvl2 = numeric(),
  ncells.lvl2 = numeric(),
  ncells.mal.lvl1 = numeric(),
  ncells.tme.lvl1 = numeric(),
  ncells.lvl1 = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(studies)) {
  seu <- readRDS(paste0(studies[i], ".Rds"))
  counts.lvl3 <- open_matrix_dir(studies[i])
  metadata.lvl3 <- seu@meta.data
  this_study <- gsub("_v5", "", basename(studies[i]))
  cell_names <- paste(colnames(counts.lvl3), i, sep = "_")
  colnames(counts.lvl3) <- cell_names
  rownames(metadata.lvl3) <- cell_names
  
  # Write level 3 BPCells matrices (raw data fater preprocessing)
  write_matrix_dir(
    mat = counts.lvl3,
    dir = paste0(
      "./v5/lvl3/",
      this_study,
      "_v5"
    )
  )
  
  data.list_lvl3[[this_study]] <- counts.lvl3
  metadata.lvl3$study <- this_study
  metadata.list_lvl3[[this_study]] <- metadata.lvl3
  
  
  # Filter cells with less than 100 malignant cells and without clinical 
  # annotation for level 2.
  metadata.lvl2 <- metadata.lvl3 %>%
    rownames_to_column("cell") %>%
    mutate("sample" = str_trim(sample)) %>%
    inner_join(y = clinical,
               by = c("sample" = "sample", "study" = "study")) %>%
    mutate(
      refined_tumor_site = case_when(refined_tumor_site == "" ~ "Unknown", 
                                     TRUE ~ refined_tumor_site)
    ) %>%
    mutate(original_cell_id = gsub(
      pattern = "\\.\\.\\..*$",
      ## annoying ... by seurat
      replacement = "",
      x = cell
    )) %>%
    column_to_rownames("original_cell_id")
  counts.lvl2 <- counts.lvl3[, rownames(metadata.lvl2)]
  
  write_matrix_dir(
    mat = counts.lvl2,
    dir = paste0(
      "./v5/lvl2/",
      this_study,
      "_v5"
    )
  )
  
  data.list_lvl2[[this_study]] <- counts.lvl2
  metadata.list_lvl2[[this_study]] <- metadata.lvl2
  
  
  # Subset of level 2 keeping only malignant cells (removing predicted tumor 
  # from brca samples)
  metadata.lvl1 <- metadata.lvl2 %>%
    filter(metadata.lvl2$malignancy == TRUE)
  counts.lvl1 <- counts.lvl2[, rownames(metadata.lvl1)]
  
  write_matrix_dir(
    mat = counts.lvl1,
    dir = paste0(
      "./v5/lvl1/",
      this_study,
      "_v5"
    )
  )
  
  data.list_lvl1[[this_study]] <- counts.lvl1
  metadata.list_lvl1[[this_study]] <- metadata.lvl1
  
  
  # Get a report of cell counts in each step.
  ncells_lvl3 <- summarize_cells(metadata.lvl3, "lvl3")
  ncells_lvl2 <- summarize_cells(metadata.lvl2, "lvl2")
  ncells_lvl1 <- summarize_cells(metadata.lvl1, "lvl1")
  
  
  ncells_all_lvls <- ncells_lvl3 %>%
    full_join(ncells_lvl2, by = "sample") %>%
    full_join(ncells_lvl1, by = "sample") %>%
    replace_na(setNames(as.list(rep(0, 9)), colnames(report_df)[3:11])) %>%
    mutate(study = this_study) %>%
    select(colnames(report_df)) %>%
    as.data.frame()
  
  report_df <- rbind(report_df, ncells_all_lvls)
  print(this_study)
}

# Save the report summary.
write.table(
  report_df,
  "./v5/report_levels.tsv",
  row.names = FALSE
)

# Merged all studies from level 3 into a single Seurat object.
metadata_lvl3 <- Reduce(rbind, metadata.list_lvl3)
seu.lvl3 <- CreateSeuratObject(counts = data.list_lvl3, meta.data = metadata_lvl3)
saveRDS(
  object = seu.lvl3,
  "./v5/lvl3/seu_lvl3.rds"
)

# Merged all studies from level 2 into a single Seurat object.
metadata_lvl2 <- Reduce(rbind, metadata.list_lvl2)
seu.lvl2 <- CreateSeuratObject(counts = data.list_lvl2, meta.data = metadata_lvl2)
saveRDS(
  object = seu.lvl2,
  "./v5/lvl2/seu_lvl2.rds"
)

# Merged all studies from level 2 into a single Seurat object.
metadata_lvl1 <- Reduce(rbind, metadata.list_lvl1)
seu.lvl1 <- CreateSeuratObject(counts = data.list_lvl1, meta.data = metadata_lvl1)
saveRDS(
  object = seu.lvl1,
  "./v5/lvl1/seu_lvl1.rds"
)

# Repair the harmonized cell annotation of level2 done by Oscar. Change 'Stem' to 
# 'Stromal' in the cell_type_main column for MSCs in the cell_type_fine column.
tcca_annot <- read.table(
  "./tcca_old/tcca_annotation_raw.tsv",
  header = TRUE,
  sep = "\t"
)

tcca_annot[tcca_annot$cell_type_fine == "MSC", "cell_type_main"] <- "Stromal"
tcca_annot <- tcca_annot[tcca_annot$cell %in% colnames(seu.lvl2), ]
tcca_annot["cell_id"] = 1:nrow(tcca_annot)
write.table(tcca_annot,
            file = "./tcca/tcca_annotation_raw.tsv",
            row.names = FALSE,
            sep = "\t")
seu.lvl2 <- readRDS(
  "./v5/lvl2/seu_lvl2.rds"
)


# Write h5ad object from BPCells joined matrix to perform later integration.
joined.lvl2 <- JoinLayers(seu.lvl2)
write_matrix_anndata_hdf5(mat = joined.lvl2[["RNA"]]$counts,
                          path = "./tcca/tcca_raw_mat.h5ad")
