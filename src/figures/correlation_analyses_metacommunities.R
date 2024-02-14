library(BPCells)
library(ComplexHeatmap)
library(tidyverse)

## Source TCCA palette
source(file = "src/figures/TCCA_palette.R")


## Correlation analysis
## Load pancancer bc mat

mat <- BPCells::open_matrix_dir(
  dir = "results/beyondcell_bp/sketch_mat_beyondcell/"
  )

## no ctrp drugs
no_ctrp <- data.table::fread("reference/final_moas - Collapsed.tsv") %>%
  filter(studies != "CTRP") %>%
  select(IDs) %>%
  distinct() %>%
  pull(IDs)

mat <- t(mat[no_ctrp, ])

correlation_mat <- cor(x = as.matrix(mat))

## load metacoms new
metacom_new <- read.table(
  "results/modules/annotated/metagroup_patients_untreated_consensus_drugs.tsv"
  ) %>%
  group_by(signature) %>%
  arrange(desc(n.appearances)) %>%
  distinct(signature, .keep_all = TRUE) %>%
  as.data.frame()

rownames(metacom_new) <- metacom_new$signature
metacom_new <- metacom_new[, c("signature", "meta_community", "collapsed.MoAs")]


## get drug names
drugs <- data.table::fread("reference/final_moas - Collapsed.tsv") %>%
  select(IDs, preferred.drug.names, collapsed.MoAs) %>%
  filter(IDs %in% no_ctrp) %>%
  mutate(
    collapsed.MoAs = case_when(
      collapsed.MoAs %in% names(MoAs_colors) ~ collapsed.MoAs,
      TRUE ~ "Other"
    )
  ) %>%
  distinct() %>%
  as.data.frame()


metacom_colors <- ecs_colors[1:5]
names(metacom_colors) <- c("1","2","3","4", "None")

moa_pals <- list(
  "Mechanism of action" = MoAs_colors,
  "Metacommunity" = metacom_colors
)

metacom_factor <- metacom_new[drugs$IDs, "meta_community"]
metacom_factor[is.na(metacom_factor)] <- "None"
names(metacom_factor) <- drugs$IDs
metacom_factor <- as.factor(metacom_factor)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  "Mechanism of action" = drugs$collapsed.MoAs,
  "Metacommunity" = metacom_factor,
  which = "column",
  col = moa_pals,
  show_annotation_name = FALSE
)

cmat <- ComplexHeatmap::Heatmap(
  matrix = correlation_mat,
  cluster_columns = FALSE,
  cluster_column_slices = TRUE,
  column_order = names(metacom_factor[order(metacom_factor)]),
  column_split = metacom_factor,
  cluster_rows = FALSE,
  cluster_row_slices = TRUE,
  row_order = names(metacom_factor[order(metacom_factor)]),
  row_split = metacom_factor,
  top_annotation = top_annotation,
  column_names_gp = gpar(fontsize = 4),
  row_names_gp = gpar(fontsize = 4)
)


## cmat with CTRP
mat <- BPCells::open_matrix_dir(
  dir = "results/beyondcell_bp/sketch_mat_beyondcell_5k"
)

mat <- t(mat)

correlation_mat <- cor(x = as.matrix(mat))

## load metacoms new
metacom_new <- read.table(
  "results/modules_full/annotated/metagroup_patients_untreated_consensus_drugs.tsv"
) %>%
  group_by(signature) %>%
  arrange(desc(n.appearances)) %>%
  distinct(signature, .keep_all = TRUE) %>%
  as.data.frame()

rownames(metacom_new) <- metacom_new$signature
metacom_new <- metacom_new[, c("signature", "meta_community", "collapsed.MoAs")]


## get drug names
drugs <- data.table::fread("reference/final_moas - Collapsed.tsv") %>%
  select(IDs, preferred.drug.names, collapsed.MoAs) %>%
  mutate(
    collapsed.MoAs = case_when(
      collapsed.MoAs %in% names(MoAs_colors) ~ collapsed.MoAs,
      TRUE ~ "Other"
    )
  ) %>%
  distinct() %>%
  as.data.frame()


metacom_colors <- ecs_colors[1:7]
names(metacom_colors) <- c("1","2","3","4", "5", "6", "None")

moa_pals <- list(
  "Mechanism of action" = MoAs_colors,
  "Metacommunity" = metacom_colors
)

metacom_factor <- metacom_new[drugs$IDs, "meta_community"]
metacom_factor[is.na(metacom_factor)] <- "None"
names(metacom_factor) <- drugs$IDs
metacom_factor <- as.factor(metacom_factor)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  "Mechanism of action" = drugs$collapsed.MoAs,
  "Metacommunity" = metacom_factor,
  which = "column",
  col = moa_pals,
  show_annotation_name = FALSE
)

cmat <- ComplexHeatmap::Heatmap(
  matrix = correlation_mat,
  cluster_columns = TRUE,
  cluster_column_slices = TRUE,
 # column_order = names(metacom_factor[order(metacom_factor)]),
 # column_split = metacom_factor,
  cluster_rows = TRUE,
 clustering_distance_rows = "pearson",
 clustering_distance_columns = "pearson",
  cluster_row_slices = TRUE,
#  row_order = names(metacom_factor[order(metacom_factor)]),
#  row_split = metacom_factor,
  top_annotation = top_annotation,
  column_split = 6,
  row_split = 6,
  column_names_gp = gpar(fontsize = 4),
  row_names_gp = gpar(fontsize = 4)
)

cmat_raw <- ComplexHeatmap::Heatmap(
  matrix = correlation_mat,
  show_row_names = FALSE,
  show_column_names = FALSE,
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  column_split = 6,
  row_split = 6
)


## repeat in mt1
clinical <- data.table::fread("results/annotation/beyondcell_with_therapeutic_clusters.tsv")

## mt1 
## Use right join to get rid of duplicate samples
clinical_module_annot <- clinical %>%
  mutate(
    metagroup = case_when(
      treated == FALSE & study != "cell_lines_gabriella_kinker" & sample_type == "p" ~ "patient_primary_untreated",
      treated == FALSE & study != "cell_lines_gabriella_kinker" & sample_type == "m" ~ "patient_metastatic_untreated",
      treated == TRUE & study != "cell_lines_gabriella_kinker" & sample_type == "p" ~ "patient_primary_treated",
      treated == TRUE & study != "cell_lines_gabriella_kinker" & sample_type == "m" ~ "patient_metastatic_treated",
      study == "cell_lines_gabriella_kinker" ~ "cell_line",
      TRUE ~ "other"
    )
  )

mat_mt1 <- BPCells::open_matrix_dir(
  dir = "results/beyondcell_bp/sketch_mat_beyondcell_5k"
)

cells_mt1 <- clinical_module_annot %>%
  filter(metagroup == "patient_primary_untreated")

cells_to_keep <- colnames(mat_mt1) %in%  cells_mt1$new_cell_id 

mat_mt1 <- mat_mt1[, cells_to_keep]

correlation_mat <- cor(x = t(as.matrix(mat_mt1)))

## repeat metacoms mt1 with ctrp using cells from mt1
drugs <- data.table::fread("reference/final_moas - Collapsed.tsv") %>%
  select(IDs, preferred.drug.names, collapsed.MoAs, studies) %>%
  mutate(
    collapsed.MoAs = case_when(
      collapsed.MoAs %in% names(MoAs_colors) ~ collapsed.MoAs,
      TRUE ~ "Other"
    )
  ) %>%
  distinct() %>%
  as.data.frame()

metacom_colors <- ecs_colors[1:7]
names(metacom_colors) <- c("1","2","3","4", "5", "6", "None")

moa_pals <- list(
    "Mechanism of action" = MoAs_colors,
    "Metacommunity" = metacom_colors
)

metacom_factor <- metacom_new[drugs$IDs, "meta_community"]
metacom_factor[is.na(metacom_factor)] <- "None"
names(metacom_factor) <- drugs$IDs
metacom_factor <- as.factor(metacom_factor)



## get coll annotation
top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  "Mechanism of action" = drugs$collapsed.MoAs,
  "Collection" = drugs$studies,
  "Metacommunity" = metacom_colors,
  which = "column",
  col = moa_pals,
  show_annotation_name = FALSE
)
cmat <- ComplexHeatmap::Heatmap(
  matrix = correlation_mat,
  cluster_columns = TRUE,
  cluster_column_slices = TRUE,
  #column_order = names(metacom_factor[order(metacom_factor)]),
  #column_split = metacom_factor,
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  #row_order = names(metacom_factor[order(metacom_factor)]),
  #row_split = metacom_factor,
  top_annotation = top_annotation,
  column_names_gp = gpar(fontsize = 4),
  row_names_gp = gpar(fontsize = 4),
  clustering_distance_columns = "pearson",
  clustering_distance_rows = "pearson"
)

cmat_raw <- ComplexHeatmap::Heatmap(
  matrix = correlation_mat,
  top_annotation = top_annotation,
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  show_row_names = FALSE,
  show_column_names = FALSE,
  clustering_distance_rows = "pearson",
  clustering_distance_columns = "pearson",
  column_split = 6,
  row_split = 6
)
ht <- draw(cmat_raw)
orders <- row_order(ht)

orders_dt_1 <- rownames(correlation_mat)[orders[[1]]]

orders[[1]] <- rownames(correlation_mat)[orders[[1]]]
orders[[2]] <- rownames(correlation_mat)[orders[[2]]]
orders[[3]] <- rownames(correlation_mat)[orders[[3]]]
orders[[4]] <- rownames(correlation_mat)[orders[[4]]]
orders[[5]] <- rownames(correlation_mat)[orders[[5]]]
orders[[6]] <- rownames(correlation_mat)[orders[[6]]]

saveRDS(object = orders, file = "metacoms_hc.rds")

