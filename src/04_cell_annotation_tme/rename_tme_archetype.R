library(dplyr)
library(tidyverse)
setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/v5")

tcca_metadata <- read.table("tcca_metadata.tsv", sep = "\t", header = TRUE)

tcca_metadata <- tcca_metadata %>%
  mutate(tme_archetype = case_when(
    tme_archetype == "Immune_rich_cDC1_bias" ~ "Immune_rich",
    tme_archetype == "Immune_rich_Treg_cDC1_cDC2_bias" ~ "Immune_rich_Treg_cDC2_bias",
    tme_archetype == "Tcell_centric_Mo_cDC1_cDC2_bias" ~ "Tcell_centric",
    tme_archetype == "Myeloid_centric" ~ "Myeloid_centric",
    tme_archetype == "Myeloid_centric_Mp_bias" ~ "Myeloid_centric_Mp_bias",
    tme_archetype == "Myeloid_centric_Mo_cDC1_cDC2_bias" ~ "Myeloid_centric_Mo_bias(high)",
    tme_archetype == "Myeloid_centric_Endo_Mo_cDC2_bias" ~ "Myeloid_centric_Mo_bias(medium)",
    tme_archetype == "Immune_stromal_cDC1_bias" ~ "Immune_stromal",
    tme_archetype == "Immune_stromal_Endolike_Treg_cDC1_cDC2_bias" ~ "Immune_stromal_Endolike",
    tme_archetype == "Immune_stromal_CAFlike_cDC1_cDC2_bias" ~ "Immune_stromal_CAFlike_Mp_bias",
    tme_archetype == "Immune_desert_CAFlike" ~ "Immune_desert_CAFlike",
    tme_archetype == "Immune_desert_cDC1_bias" ~ "Immune_stromal_desert",
    TRUE ~ tme_archetype
  ),
  tme_split = str_split_fixed(tme_archetype, "_", n = 3),
  tme_archetype_group = case_when(
      tme_archetype == "Immune_stromal_desert" ~ "Immune_desert",
      tme_split[,3] == "" ~ tme_archetype,
      TRUE ~ paste(tme_split[,1], tme_split[,2], sep = "_")
    )
  ) %>%
  select(-tme_split) %>%
  relocate(tme_archetype_group, .before = tme_archetype)

write.table(tcca_metadata, "tcca_metadata.tsv", sep = "\t")
