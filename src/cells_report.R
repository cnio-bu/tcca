library("tidyverse")
library("patchwork")

all_reports <- list.files(
    path = "results/non_malignant_assessing_test/reports",
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
            grepl(pattern = "bone", x = study) ~ "Bone sarcoma",
            grepl(pattern = "breast", x = study) ~ "Breast cancer",
            grepl(pattern = "brmets", x = study) ~ "Brain metastasis",
            grepl(pattern = "cc", x = study) ~ "Cervical cancer",
            grepl(pattern = "cll", x = study) ~ "Chronic lymphocytic leukemia",
            grepl(pattern = "gbm", x = study) ~ "Glioma/Glioblastoma",
            grepl(pattern = "luad", x = study) ~ "Lung adenocarcinoma",
            grepl(pattern = "pancancer", x = study) ~ "Pancancer",
            grepl(pattern = "pdac", x = study) ~ "Pancreatic adenocarcinoma",
            grepl(pattern = "pleural", x = study) ~ "Pleuroblastoma",
            grepl(pattern = "rcell", x = study) ~ "Renal cell carcinoma",
            grepl(pattern = "synovial", x = study) ~ "Synovial sarcoma",
            grepl(pattern = "urothelial", x = study) ~ "Urothelial carcinoma",
            grepl(pattern = "mmieloma", x = study) ~ "Multiple mieloma"
            
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
        
        
        
report_summary_arranged <- full_report_annotated %>%
    group_by(cancer_type) %>%
    mutate(
        avg_cells_by_type = mean(malignants)
    ) %>%
    arrange(desc(avg_cells_by_type))
   
report_summary_arranged$cancer_type <- fct_reorder(
    report_summary_arranged$cancer_type,
    report_summary_arranged$avg_cells_by_type
    )

## plot averages
avg_malignants_by_tumor <- ggplot(data = report_summary_arranged) + 
    geom_boxplot(aes(y=cancer_type, x = malignants)) +
    scale_x_continuous(n.breaks = 10) +
    xlab(label = "Malignant cells") +
    ylab(label = "") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    ) +
    ggtitle(label = "Average proportion of malignant cells by sample and tumor type")


report_summary_all_malignants <- full_report_annotated %>%
    group_by(cancer_type) %>%
    summarise(
        n.malignants = sum(malignants),
        n.cells = sum(cells),
        n.samples = n()
    )  %>%
    arrange(n.malignants) %>%
    mutate(
        cancer_type = as_factor(cancer_type)
    )

## plot totals
malignants_by_tumor <- ggplot(data=report_summary_all_malignants, aes(x = n.malignants, y = cancer_type)) + 
    geom_bar(stat = "identity") +
    geom_text(stat="identity", aes(label = n.malignants), hjust = -1) +
    scale_x_continuous(n.breaks = 10) +
    xlab(label = "Malignant cells") +
    ylab(label = "") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    ) +
    ggtitle(label = "Total proportion of malignant cells by tumor type")

## plots by proportion
report_summary_all_malignants <- report_summary_all_malignants %>%
    mutate(
        prop.total = n.malignants / sum(n.malignants),
        prop.total = round(prop.total, digits = 3)
    )


## plot props
malignants_by_tumor <- ggplot(data=report_summary_all_malignants, aes(x = prop.total, y = cancer_type)) + 
    geom_bar(stat = "identity") +
    scale_x_continuous(n.breaks = 18, labels = scales::percent_format()) +
    xlab(label = "Malignant cells") +
    ylab(label = "") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    ) +
    ggtitle(label = "Total proportion of malignant cells by tumor type")


report_summary_all_malignants <- full_report_annotated %>%
    group_by(cancer_type) %>%
    summarise(
        n.malignants = sum(malignants),
        n.cells = sum(cells),
        n.samples = n()
    )  %>%
    arrange(n.samples) %>%
    mutate(
        cancer_type = as_factor(cancer_type)
    )

## samples by tumor
samples_by_tumor <- ggplot(data = report_summary_all_malignants, aes(x = n.samples, y = cancer_type)) +
    geom_bar(stat = "identity") +
    scale_x_continuous(n.breaks = 10) +
    xlab(label = "Samples") +
    ylab(label = "") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    ) +
    ggtitle(label = "Samples by tumor type")


ggsave(
    plot = samples_by_tumor,
    filename = "results/samples_by_tumor_type.png",
    dpi = 300,
    height = 9,
    width = 16
    )


ggsave(
    plot = malignants_by_tumor,
    filename = "results/proportion_of_malignants_per_tumor.png",
    dpi = 300,
    height = 9,
    width = 16
)

study_wise_report <- full_report_annotated %>%
    group_by(cancer_type) %>%
    summarise(
        n.malignants = sum(malignants),
        n.cells = sum(cells),
        n.samples = n(),
        n.studies = sum(length(unique(study)))
    )  %>%
    arrange(n.studies) %>%
    mutate(
        cancer_type = as_factor(cancer_type)
    ) 

## studies by tumor
samples_by_tumor <- ggplot(data = study_wise_report, aes(x = n.studies, y = cancer_type)) +
    geom_bar(stat = "identity") +
    scale_x_continuous(n.breaks = 4) +
    xlab(label = "Samples") +
    ylab(label = "") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    ) +
    ggtitle(label = "Studies by tumor type")


ggsave(
    plot = samples_by_tumor,
    filename = "results/studies_per_tumor.png",
    dpi = 300,
    height = 9,
    width = 16
)
