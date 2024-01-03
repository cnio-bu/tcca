library(cluster)
library(data.table)
library("tidyverse")
library("ComplexHeatmap")

clinical_database <- data.table::fread(
  "results/annotation/clinical_metadata_v4_clean.tsv"
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
modules <- lapply(
    module_list,
    data.table::fread, select = c(
        "bicluster",
        "cluster_contribution",
        "signature",
        "information_content"
        )
    )

modules <- modules %>%
  map(distinct) %>%
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

write_tsv(
    x = modules_annotated_clinical,
    file = "results/annotation/modules_annotated_clinical.tsv"
    )

