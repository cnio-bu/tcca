library("tidyverse")
library("patchwork")

all_reports <- list.files(
    path = "results/reports/",
    full.names = TRUE
    )

full_report <- all_reports %>%
    map(read_tsv, id = "study") %>%
    bind_rows() 


full_report_annotated <- full_report %>%
    mutate(
        study = basename(study),
        study = stringr::str_remove_all(string = study, pattern = "cells_"),
        study = stringr::str_remove_all(string = study, pattern = ".tsv"),
    ) %>%
    mutate(
        cancer_type = case_when(
            grepl(pattern = "adrenalnb", x = study) ~ "Adrenal Neuroblastoma",
            grepl(pattern = "aml", x = study) ~ "Acute Myeloid Leukemia",
            grepl(pattern = "all", x = study) ~ "Acute lymphocytic Leukemia",
            grepl(pattern = "bone", x = study) ~ "Bone sarcoma",
            grepl(pattern = "breast", x = study) ~ "Breast cancer",
            grepl(pattern = "brca", x = study) ~ "Breast cancer",
            grepl(pattern = "brmets", x = study) ~ "Brain metastasis",
            grepl(pattern = "cc", x = study) ~ "Cervical cancer",
            grepl(pattern = "cll", x = study) ~ "Chronic lymphocytic leukemia",
            grepl(pattern = "crc", x = study) ~ "Colon adenocarcinoma",
            grepl(pattern = "esca", x = study) ~ "Esophageal cancer",
            grepl(pattern = "gbm", x = study) ~ "Glioma/Glioblastoma",
            grepl(pattern = "luad", x = study) ~ "Lung adenocarcinoma",
            grepl(pattern = "pancancer", x = study) ~ "Pancancer",
            grepl(pattern = "pdac", x = study) ~ "Pancreatic adenocarcinoma",
            grepl(pattern = "pleural", x = study) ~ "Pleuroblastoma",
            grepl(pattern = "rcell", x = study) ~ "Renal cell carcinoma",
            grepl(pattern = "synovial", x = study) ~ "Synovial sarcoma",
            grepl(pattern = "urothelial", x = study) ~ "Urothelial carcinoma",
            grepl(pattern = "mmieloma", x = study) ~ "Multiple mieloma",
            grepl(pattern = "uvm", x = study) ~ "Uveal melanoma",
            grepl(pattern = "nsclc", x = study) ~ "Non small cell lung cancer",
            grepl(pattern = "prad", x = study) ~ "Prostate adenocarcinoma",
            grepl(pattern = "skcm", x = study) ~ "Skin cutaneous melanoma",
            grepl(pattern = "cell_lines", x = study) ~ "Cancer cell line"
            
        )
    ) %>%
    mutate(
        is_pancancer = cancer_type == "Pancancer",
        cancer_type = case_when(
            cancer_type == "Pancancer" & grepl(pattern = "OV", x = sample) ~ "Ovarian cancer",
            cancer_type == "Pancancer" & grepl(pattern = "BRCA", x = sample) ~ "Breast cancer",
            cancer_type == "Pancancer" & grepl(pattern = "PDAC", x = sample) ~ "Pancreatic adenocarcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "COAD", x = sample) ~ "Colon adenocarcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "LIHC", x = sample) ~ "Liver hepatocellular carcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "CHCA", x = sample) ~ "Cholangiocarcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "PRAD", x = sample) ~ "Prostate adenocarcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "KIRCH", x = sample) ~ "Renal cell carcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "UCECH", x = sample) ~ "Endometrial carcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "LUAD", x = sample) ~ "Lung adenocarcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "HNSCH", x = sample) ~ "Head and neck squamous",
            cancer_type == "Pancancer" & grepl(pattern = "SKSCH", x = sample) ~ "Skin squamous cell carcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "GISTH", x = sample) ~ "Gastrointestinal stromal cancer",
            cancer_type == "Pancancer" & grepl(pattern = "GBM", x = sample) ~ "Glioblastoma",
            cancer_type == "Pancancer" & grepl(pattern = "OGDH", x = sample) ~ "Oligodendroglioma",
            cancer_type == "Pancancer" & grepl(pattern = "GBM", x = sample) ~ "Glioblastoma",
            cancer_type == "Pancancer" & grepl(pattern = "GBM", x = sample) ~ "Glioblastoma",
            cancer_type == "Pancancer" & grepl(pattern = "BC-", x = sample) ~ "Breast cancer",
            cancer_type == "Pancancer" & grepl(pattern = "PC-", x = sample) ~ "Prostate adenocarcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "M-", x = sample) ~ "Lymph. metastasis",
            cancer_type == "Pancancer" & grepl(pattern = "sc5rJUQ", x = sample) ~ "Breast cancer",
            cancer_type == "Pancancer" & grepl(pattern = "scrEXT", x = sample) ~ "Colon adenocarcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "BT", x = sample) ~ "Lung adenocarcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "scrBT", x = sample) ~ "Lung adenocarcinoma",
            cancer_type == "Pancancer" & grepl(pattern = "scrSOL", x = sample) ~ "Ovarian cancer",
            TRUE ~ cancer_type
        )
    ) %>%
    mutate(
        tumor_type = case_when(
            str_detect(string = cancer_type, pattern = "metastasis") ~ "Metastatic",
            TRUE ~ "Primary"
        )
    )
        


## keep those with > 100 malignants
samples_to_keep <- full_report_annotated %>%
    filter(malignants >= 100)

cells_by_study <- samples_to_keep %>%
    group_by(study) %>%
    summarise(
        n.samples = n(),
        n.cells = sum(malignants)
    )


samples_by_study <- samples_to_keep %>%
    group_by(study) %>%
    summarise(
        sample = sample
    )

write.csv(x = samples_by_study, file = "results/all_samples_annotated.tsv")
