library(cluster)
library(data.table)
library("tidyverse")
library("ComplexHeatmap")

clinical_database <- data.table::fread(
  "results/annotation/clinical_metadata_v2_clean.tsv"
) %>%
  filter(
    tumor_subtype != "predicted_tumour" ## Get rid of non clear tumoral stuff
  )

module_list <- list.files(
  "results/modules",
  full.names = TRUE,
  recursive = TRUE,
  pattern = ".tsv"
  )

## Load all modules
modules <- module_list %>%
  map(data.table::fread) %>%
  bind_rows(.id = "sample_study")

modules$sample_study <- module_list[as.numeric(modules$sample_study)]

modules_annotated <- modules %>%
  mutate(
    sample_study = stringr::str_remove(
      pattern = "results/modules/", string = sample_study)
  ) %>%
  separate(
    col = sample_study,
    into = c("study", "sample"),
    sep = "/"
  ) %>%
  mutate(
    sample = stringr::str_remove(string = sample, pattern = "_clusters.tsv")
  )

## Use right join to get rid of duplicate samples
modules_annotated_clinical <- modules_annotated %>%
  right_join(
    y = clinical_database,
    by = c("sample" = "sample", "study" = "study")
  ) %>%
  mutate(
    metagroup = case_when(
      treated == "f" & study != "cell_lines_gabriella_kinker" & sample_type == "p" ~ "patient_primary_untreated",
      treated == "f" & study != "cell_lines_gabriella_kinker" & sample_type == "m" ~ "patient_metastatic_untreated",
      treated == "t" & study != "cell_lines_gabriella_kinker" & sample_type == "p" ~ "patient_primary_treated",
      treated == "t" & study != "cell_lines_gabriella_kinker" & sample_type == "m" ~ "patient_metastatic_treated",
      study == "cell_lines_gabriella_kinker" ~ "cell_line",
      TRUE ~ "other"
    )
  )

## generate a summary for F/G
module_summary <- modules_annotated_clinical %>%
  group_by(metagroup) %>%
  summarise(
    n.samples = length(unique(sample)),
    n.tumor_types = length(unique(tumor_type))
  )


## Perform aggregation of metagroup 1
modules_mt1 <- modules_annotated_clinical %>%
  filter(metagroup == "patient_primary_untreated")

## 29 tumor types in MT1 ppu
module_mt1_summary <- modules_mt1 %>%
  group_by(tumor_type) %>%
  summarise(
    n.samples = length(unique(sample))
  )

## Keep tumor types with >= 5 samples
enough_samples <- module_mt1_summary %>%
  filter(n.samples >= 5) %>%
  pull(tumor_type)


## tumor level aggregation
for(tumor in enough_samples){
  this_subset <- modules_mt1 %>%
    filter(
      tumor_type == tumor
    )
  
  ## Generate module comparison matrix
  
}

test <- modules_mt1 %>%
    filter(tumor_type == "BRCA")

test_wide_2 <- test %>%
    select(sample, signature, bicluster) %>%
    filter(bicluster <= 5) %>%
    mutate(
        sample = paste(sample, bicluster, sep = "_")
    )  %>%
    group_by(sample, signature) %>%
    tally() %>%
    spread(signature, n, fill = 0) %>%
    as.data.frame()

rownames(test_wide_2) <- test_wide_2$sample
test_wide_2$sample <- NULL

test_wide_2 <- as.matrix(test_wide_2)
test_wide_2 <- ifelse(test_wide_2 == 0,FALSE,TRUE)

test4 <- proxy::dist(
    x = test_wide_2,
    by_rows = FALSE,
    pairwise = TRUE,
    method = "jaccard",
    )

clust_mat <- hclust(test4, method = "complete")

library(beyondcell)

b <- drugInfo$IDs
b <- b %>%
    filter(collections == "SSc")

moas <- drugInfo$MoAs
b_annot <- b %>%
    left_join(y = moas, by = "IDs")

egfrs <- b_annot %>%
    filter(main.MoAs %in% c("Cell cycle arrest","DNA replication inhibitor")) %>%
    pull("IDs")


egfrs_inhit <- data.frame("drug" = names(test4), "MoA" = names(test4) %in% egfrs)
egfrs_inhit$MoA <- as.factor(egfrs_inhit$MoA)
levels(egfrs_inhit$MoA) <- c("Other", "Cell cycle arrest")

rownames(egfrs_inhit) <- egfrs_inhit$drug
egfrs_inhit$drug <- NULL

egfr_annot <- HeatmapAnnotation(
    df = egfrs_inhit, name = "DNA replication inhibitor", which = "column"
)


heat <- Heatmap(
    matrix = 1 - as.matrix(test4),
    cluster_rows = clust_mat,
    cluster_columns = clust_mat,
    show_column_names = TRUE,
    show_column_dend = FALSE,
    show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 2),
    top_annotation = egfr_annot
    )


