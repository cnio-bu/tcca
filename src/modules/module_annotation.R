library("tidyverse")

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
      treated == "t" & study != "cell_lines_gabriella_kinker" & sample_type == "p" ~ "patient_metastatic_treated",
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
combinations <- apply(combn(unique(test$sample), 2), 2, function(z) paste(sort(z), collapse = "-"))

test2 <- test %>% 
  group_by(bicluster, signature) %>%
  summarise(
    names = paste(sort(unique(sample)), collapse = "-")
    ) %>%
  right_join(tibble(names = combinations), by = "names") %>%
  group_by(names, bicluster) %>%
  summarise(
    n_overlap = sum(!is.na(signature))
    )
