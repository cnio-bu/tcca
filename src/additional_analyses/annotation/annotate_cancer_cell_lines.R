library("tidyverse")

cell_line_annotation <- data.table::fread("raw/cell_line_annotation.csv")

cancer_lines_profiled <- data.table::fread("results/all_malignant_samples.tsv") %>%
  select(-V1) %>%
  filter(
    study == "cell_lines_gabriella_kinker",
    malignants >= 100
  ) %>%
  select(sample)


cancer_lines_profiled_annotated <- cancer_lines_profiled %>%
  left_join(
    y = cell_line_annotation,
    by = c("sample" = "CCLE_Name")
  ) %>%
  select(
    sample,
    sex,
    age,
    primary_disease,
    primary_or_metastasis,
    lineage,
    lineage_subtype,
    lineage_sub_subtype,
    lineage_molecular_subtype,
    sample_collection_site
  ) %>%
  arrange(sample)


write.csv(
  x = cancer_lines_profiled_annotated,
  file = "results/cell_lines_annotation/cell_lines_profiled_database.csv"
  )
