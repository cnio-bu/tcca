library(dplyr)
library(tidyverse)

setwd("/Users/mariagb/Library/CloudStorage/OneDrive-CNIO/2nd_year/bc-meta/")

# Concordance of approved drugs for cancer type
drugs <- read.table("code/Shiny_TCCA_App_final/www/drug_response_subclone_final.tsv", sep = "\t", header = TRUE)

# Adjutsted dictionary with exact values from Approved.Cancer.Types
tumor_type_map <- list(
  "ALL"  = c("acute lymphoblastic leukemia", "lymphoblastic leukemia"),
  "LAML" = c("acute myeloid leukemia", "acute nonlymphocytic leukemia"),
  "BCC"  = c("basal cell carcinoma"),
  "SARC" = c("sarcoma", "gastrointestinal stromal tumor"),
  "BRCA" = c("breast cancer"),
  "SKCM" = c("melanoma"),
  "SCLC" = c("small cell lung cancer"),
  "LUAD" = c("non-small cell lung cancer"),
  "LUSC" = c("non-small cell lung cancer"),
  "LCLC" = c("non-small cell lung cancer"),
  "OV"   = c("ovarian cancer", "fallopian tube cancer", "peritoneal cancer"),
  "COAD" = c("colorectal cancer", "adenocarcinoma of the colon and rectum",
              "carcinoma of the colon or rectum"),
  "READ" = c("colorectal cancer", "adenocarcinoma of the colon and rectum",
              "carcinoma of the colon or rectum"),
  "BLCA" = c("carcinoma in situ of the urinary bladder"),
  "CESC" = c("cervical cancer"),
  "CHOL" = c("cholangiocarcinoma"),
  "CLL"  = c("chronic lymphocytic leukemia", "small lymphocytic lymphoma"),
  "ESCA" = c("adenocarcinoma of the gastroesophageal junction",
              "gastric adenocarcinoma", "esophageal carcinoma"),
  "GBM" = c("astrocytoma"), 
  "LGG" = c("astrocytoma"),
  "MM" = c("multiple myeloma", "mantle cell lymphoma"),
  "UCEC" = c("endometrial carcinoma"),
  "PAAD" = c("adenocarcinoma of the pancreas", "pancreatic adenocarcinoma",
              "progressive neuroendocrine tumors of pancreatic origin"),
  "LIHC" = c("hepatocellular carcinoma"),
  "PRAD" = c("prostate cancer"),
  "KIRC" = c("renal cell carcinoma", "renal-cell carcinoma"),
  "HNSC" = c("squamous cell carcinoma of the head and neck"),
  "SKSC" = c("Skin Squamous Cell Carcinoma"), # there is not specific indication for SKSC in PanDrugs
  "UVM" = c("Uveal Melanoma") # uveal melanoma does not have approved clinical indications
)

# Matching function
is_concordant <- function(tumor_type, approved_cancers) {
  if (is.na(approved_cancers) || approved_cancers == "") return(FALSE)
  approved_list <- tolower(trimws(strsplit(approved_cancers, ",")[[1]]))
  synonyms <- tolower(tumor_type_map[[tumor_type]])
  if (is.null(synonyms)) return(FALSE)
  any(sapply(synonyms, function(s) any(grepl(s, approved_list, fixed = TRUE))))
}

# Universe: all drug-cancer type possible combinations within the 158 drugs
universe <- drugs %>%
  filter(Drug.Status == "Approved") %>%
  distinct(Drug.Name, Approved.Cancer.Types) %>%
  # For each drug, expand to the cancer types it is approved
  rowwise() %>%
  mutate(approved_for = list(tolower(trimws(strsplit(Approved.Cancer.Types, ",")[[1]])))) %>%
  ungroup()

# For each cancer type, how many out of the 158 drugs are approved for the specific cancer type
approved_per_tumor <- data.frame()

for (tt in names(tumor_type_map)) {

  syns <- tolower(tumor_type_map[[tt]])

  if (is.null(syns) || length(syns) == 0) {
    n <- 0
  } else {

    n <- universe %>%
      rowwise() %>%
      filter(
        any(
          sapply(
            syns,
            function(s)
              any(grepl(s, approved_for, fixed = TRUE))
          )
        )
      ) %>%
      nrow()
  }

  approved_per_tumor <- rbind(
    approved_per_tumor,
    data.frame(
      Refined.Tumor.Type = tt,
      n_approved_in_universe = n
    )
  )
}

# Numerator: predicted and approved for the specific tumor type
concordance <- drugs %>%
  filter(Drug.Status == "Approved") %>%
  distinct(Refined.Tumor.Type, Refined.Tumor.Type.Name, Drug.Name, Approved.Cancer.Types) %>%
  rowwise() %>%
  mutate(concordant = is_concordant(Refined.Tumor.Type, Approved.Cancer.Types)) %>%
  ungroup() %>%
  group_by(Refined.Tumor.Type, Refined.Tumor.Type.Name) %>%
  summarise(n_concordant = sum(concordant), .groups = "drop") %>%
  left_join(approved_per_tumor, by = "Refined.Tumor.Type") %>%
  mutate(pct_concordant = ifelse(n_approved_in_universe == 0, NA,
                                  n_concordant / n_approved_in_universe * 100)) 

plot_df <- concordance %>%
  filter(n_approved_in_universe > 0) %>%
  mutate(
    label = paste0(Refined.Tumor.Type, " (n=", n_approved_in_universe, ")")
  ) %>%
  arrange(
    pct_concordant,
    n_approved_in_universe,
    Refined.Tumor.Type
  )  %>%
  mutate(
    label = factor(label, levels = unique(label))
  )

# Visualize
library(ggplot2)
library(dplyr)

plot_concordance <- ggplot(
    plot_df,
    aes(
      x = label,
      y = pct_concordant,
      fill = pct_concordant
    )
  ) +
  geom_col(width = 0.8, color = "grey30", linewidth = 0.2) +

  geom_text(
    aes(label = sprintf("%.1f%%", pct_concordant)),
    hjust = -0.15,
    size = 3.5
  ) +

   scale_fill_gradientn(
  colours = c("#74A9CF", "#F7F7F7", "#F4A582"),
  name = "Concordance (%)"
) +

  coord_flip(clip = "off") +

  scale_y_continuous(
    limits = c(0, 110),
    expand = expansion(mult = c(0, 0.02))
  ) +

  labs(
    x = NULL,
    y = "Concordant approved drug predictions (%)"
  ) +

  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "right",
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  )

ggsave(
    "figures/sctherapy/concordance_approved_drugs.pdf", 
    width = 6, 
    height = 6
    )
