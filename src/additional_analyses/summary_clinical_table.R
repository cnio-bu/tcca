library(dplyr)
library(tidyverse)
library(gtsummary)
library(gt)
setwd("/home/lmgonzalezb/Documents/bc-meta/cohort_statistics/")

metadata <- read.table("tcca_metadata.tsv", sep = "\t", header = TRUE)
clinical <- metadata %>%
  select(3, 9:29, 33) %>%
  distinct()

clinical_summary <- clinical %>%
  group_by(refined_tumor_type) %>%
  summarise(
    n_patients = n(),
    n_samples = n(),
    
    female_n = sum(sex == "f", na.rm = TRUE),
    male_n = sum(sex == "m", na.rm = TRUE),
    unknown_sex_n = sum(is.na(sex) | !sex %in% c("f", "m")),
    
    primary_n = sum(sample_type == "p", na.rm = TRUE),
    metastasis_n = sum(sample_type == "m", na.rm = TRUE),
    unknown_sample_type_n = sum(is.na(sample_type) |
                                  !sample_type %in% c("p", "m")),
    
    min_age = min(age, na.rm = TRUE),
    max_age = max(age, na.rm = TRUE),
    mean_age = round(mean(age, na.rm = TRUE), 2),
    
    treated_n = sum(treated == "t", na.rm = TRUE),
    untreated_n = sum(treated == "f", na.rm = TRUE),
    unknown_treatment_n = sum(is.na(treated) |
                                !treated %in% c("t", "f"))
  ) %>%
  mutate(
    age_range = case_when(
      is.infinite(min_age) | is.infinite(max_age) ~ "—",
      min_age == max_age ~ as.character(min_age),
      TRUE ~ paste0(mean_age, " (", min_age, " - ", max_age, ")")
    ),
    sex = paste0(female_n, " / ", male_n, " / ", unknown_sex_n),
    sample_type = paste0(primary_n, " / ", metastasis_n, " / ", unknown_sample_type_n),
    treatment = paste0(treated_n, " / ", untreated_n, " / ", unknown_treatment_n)
  ) %>%
  select(refined_tumor_type,
         n_patients,
         n_samples,
         sex,
         age_range,
         sample_type,
         treatment)

clinical_summary %>%
  gt() %>%
  cols_label(
    refined_tumor_type = "Cancer Type",
    n_patients = "Total patients", 
    n_samples = "Total samples",
    sex = "Sex (F/M/U)",
    age_range = "Age (Mean [Min - Max])",
    treatment = "Treated"
  ) %>%
  tab_header(
    title = "Clinical Summary by Cancer Type"
  )
