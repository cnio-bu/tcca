library(dplyr)
library(tidyverse)
library(circlize)
set.seed(123)
setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/")

drugs <- read.table("sctherapy/results/drug_response_subclone_sctherapy.tsv", sep = "\t", header = TRUE)

# Load metadata with clinical data with mutations
clinical <- read.table("seurat/v5/clinical_metadata_v4_clean_new.tsv", sep = "\t", header = TRUE)
clinical <- clinical %>%
    mutate(study_sample = paste0(study, "_", sample)) %>%
    select(study_sample, tumor_subtype)

drugs <- drugs %>%
    mutate(study_sample = paste0(Study, "_", Sample)) %>%
    left_join(clinical, by = "study_sample") 

# Function
library(dplyr)
library(tidyr)


# FUNCTION
calculate_drug_response_pct <- function(drugs_df,
                                        drug_list,
                                        mutant_subclones,
                                        wt_subclones,
                                        status_name = "status") {

  # Mutational status per subclone
  subclone_status <- tibble(
    Subclone.Name = c(mutant_subclones, wt_subclones),
    status = c(
      rep("mutant", length(mutant_subclones)),
      rep("wt", length(wt_subclones))
    )
  )

  # Binary table: subclone × drug
  response_matrix <- expand.grid(
    Subclone.Name = subclone_status$Subclone.Name,
    Drug.Name = drug_list,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    left_join(
      drugs_df %>%
        distinct(Subclone.Name, Drug.Name) %>%
        mutate(predicted = TRUE),
      by = c("Subclone.Name", "Drug.Name")
    ) %>%
    mutate(predicted = replace_na(predicted, FALSE)) %>%
    left_join(subclone_status, by = "Subclone.Name")

  # Summary
  response_matrix %>%
    group_by(Drug.Name, status) %>%
    summarise(
      n_predicted = sum(predicted),
      n_total = n(),
      pct_predicted = round(100 * mean(predicted), 1),
      .groups = "drop"
    ) %>%
    arrange(Drug.Name, status)
}

# Drugs
# EGFR inhibitors (direct — TKIs de EGFR/HER2)
egfr_drugs <- c("AFATINIB",      # EGFR/HER2 TKI 2ª gen irreversible
                "OSIMERTINIB",   # EGFR TKI 3ª gen, T790M y activator mutations
                "NERATINIB",     # EGFR/HER2 TKI irreversible
                "LAPATINIB",     # EGFR/HER2 TKI 1ª gen
                "PELITINIB",     # pan-EGFR TKI irreversible
                "VANDETANIB")    # EGFR/VEGFR/RET TKI

# KRAS downstream inhibitors (MEK/ERK — not direct inhibitors of KRAS in drug table)
kras_drugs <- c("TRAMETINIB",    # MEK1/2 inhibitor
                "COBIMETINIB",   # MEK1 inhibitor
                "SELUMETINIB",   # MEK1/2 inhibitor
                "PD-0325901",    # MEK1/2 inhibitor
                "TAK-733",       # MEK1/2 inhibitor
                "MEK162",        # MEK1/2 inhibitor (binimetinib)
                "RO-4987655",    # MEK1/2 inhibitor
                "AS-703026",     # MEK1/2 inhibitor
                "REFAMETINIB",   # MEK1 inhibitor
                "AZD8330",       # MEK1/2 inhibitor
                "BVD-523")       # ERK1/2 inhibitor (downstream of MEK)

# BRAF inhibitors
braf_drugs <- c("VEMURAFENIB",   # BRAF V600E inhibitor
                "DABRAFENIB",    # BRAF V600E inhibitor
                "RAF265",        # pan-RAF inhibitor
                "AZ-628")        # pan-RAF inhibitor


# EGFR mutation
egfr_subclones_mut <- drugs %>%
  filter(tumor_subtype %in% c("egfr_19deletion",
                              "egfr_mutant")) %>%
  pull(Subclone.Name) %>%
  unique()

egfr_subclones_wt <- drugs %>%
  filter(tumor_subtype == "egfr_wt") %>%
  pull(Subclone.Name) %>%
  unique()

egfr_results <- calculate_drug_response_pct(
  drugs_df = drugs,
  drug_list = egfr_drugs,
  mutant_subclones = egfr_subclones_mut,
  wt_subclones = egfr_subclones_wt
)


# KRAS mutation
unique(drugs$tumor_subtype[grepl("kras", drugs$tumor_subtype, ignore.case = TRUE)])

kras_subclones_mut <- drugs %>%
  filter(tumor_subtype %in% c(
    "kras",
    "acinar_kras_mutant",
    "solid_kras_mutant",
    "mucinous_kras_mutant",
    "lepidic_kras_mutant",
    "papillary_kras_mutant",
    "kras_mutant"
  )) %>%
  pull(Subclone.Name) %>%
  unique()

kras_subclones_wt <- drugs %>%
  filter(tumor_subtype == "kras_wt") %>%
  pull(Subclone.Name) %>%
  unique()

kras_results <- calculate_drug_response_pct(
  drugs_df = drugs,
  drug_list = kras_drugs,
  mutant_subclones = kras_subclones_mut,
  wt_subclones = kras_subclones_wt
)

# BRAF mutation (only brain metastasis samples from melanoma)
braf_subclones_mut <- drugs %>%
  filter(tumor_subtype %in% c(
    "braf_r178",
    "braf_v600e"
  )) %>%
  pull(Subclone.Name) %>%
  unique()

braf_subclones_wt <- drugs %>%
  filter(tumor_subtype == "braf_wt") %>%
  pull(Subclone.Name) %>%
  unique()

braf_results <- calculate_drug_response_pct(
  drugs_df = drugs,
  drug_list = braf_drugs,
  mutant_subclones = braf_subclones_mut,
  wt_subclones = braf_subclones_wt
)

braf_results

# Dumbell plot to summarize EGFR and KRAS results
final_res <- rbind(egfr_results, kras_results)
library(ggplot2)

final_res$Drug.Name <- factor(
  final_res$Drug.Name,
  levels = rev(unique(final_res$Drug.Name))
)


plot <- ggplot(final_res, aes(x = pct_predicted, y = Drug.Name)) +
  geom_line() +
  geom_point(aes(color = status), size = 3) +
  scale_fill_manual(
    values = c(
      mutant = "#D55E00",
      wt = "#0072B2"
    ),
    labels = c(
      mutant = "Mutant",
      wt = "WT"
    )
  ) +
  labs(
    x = "Subclones with predicted response (%)",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(
    colour = "#E6E6E6",
    linewidth = 0.4
    ),
    panel.grid.minor.x = element_line(
    colour = "#F2F2F2",
    linewidth = 0.3
    ),
    panel.grid.major.y = element_line(
    colour = "#E6E6E6",
    linewidth = 0.4
    ),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_text(
      size = 10,
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 10,
      colour = "black"
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    plot.margin = margin(
      t = 5,
      r = 15,
      b = 5,
      l = 5
    )) +
  theme(legend.position = "bottom")

ggsave("sctherapy/results/egfr_kras_mutant.pdf", plot = plot, dpi = 300, width = 8, height = 8)