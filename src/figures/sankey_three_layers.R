library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggsankey)

setwd("/Users/mariagb/Library/CloudStorage/OneDrive-CentroNacionaldeInvestigacionesOncológicas/2nd_year/bc-meta/")
source("bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")
# Load cell level annotation
tcca_metadata <- read.table("tcca_metadata.tsv", sep = "\t", header = TRUE)

# Load functional MPs more active in each subclone
subclone <- read.table("functional_mps/sample_wise/subclone_top_mps.tsv", sep = "\t", header = TRUE)

# Select variables to plot
subclone <- tcca_metadata %>%
  filter(malignancy == "True" & !(scevan_subclone %in% c("", "non_tumor"))) %>%
  select(scevan_subclone, tme_archetype_group, refined_tumor_type) %>%
  distinct() %>% 
  left_join(select(subclone, scevan_subclone, top_MP_clean), by = "scevan_subclone") %>%
  rename(functional_metaprogram = top_MP_clean) %>%
  mutate(
    # Extrae todo antes del primer punto
    functional_metaprogram_family = str_extract(functional_metaprogram, "^[^.]+"),
    # Sobrescribe según categorías generales
    functional_metaprogram_family = case_when(
      str_detect(functional_metaprogram, "LineageSpecific.Hemato") ~ "LineageSpecific.Hemato",
      str_detect(functional_metaprogram, "LineageSpecific.Neural") ~ "LineageSpecific.Neural",
      str_detect(functional_metaprogram, "Secretory|Melanocyte") ~ "LineageSpecific.Other",
      TRUE ~ functional_metaprogram_family
    )
  )




# Count subclones for each combination
df_counts <- subclone %>%
  group_by(tme_archetype_group, refined_tumor_type, functional_metaprogram_family) %>%
  summarise(value = n(), .groups = "drop")

# Transform to long format
df_long <- df_counts %>%
  make_long(tme_archetype_group, refined_tumor_type, functional_metaprogram_family, value = "value")
df_long <- df_long %>%
  mutate(
    x = factor(x, levels = c("tme_archetype_group", "refined_tumor_type", "functional_metaprogram_family")),
    node = factor(node, levels = unique(c(node, next_node)))
  )
all_nodes <- unique(c(df_long$node, df_long$next_node))

# Revisar si faltan colores
missing_colors <- setdiff(all_nodes, names(c(tme_group_colors, tumor_type_colors, mp_family_colors)))
missing_colors  # si no está vacío, asigna un color


# Sankey plot
sankey <- ggplot(df_long, aes(x = x, next_x = next_x, node = node, next_node = next_node, value = value)) +
                geom_sankey(aes(fill = node), flow.alpha = 0.8) + 
                geom_sankey_label(aes(label = node), size = 3, color = "white") +
                scale_fill_manual(values = c(tme_group_colors, tumor_type_colors, mp_family_colors)) + 
                scale_x_discrete(limits = c("tme_archetype_group", "refined_tumor_type", "functional_metaprogram_family")) +
                theme_sankey(base_size = 12) + 
                geom_sankey_label(aes(label = node), size = 3, color = "white", fill = NA) +
                theme(legend.position = "none")


ggsave("sankey_plot_three_layers.png", plot = sankey, width = 15, height = 10)
