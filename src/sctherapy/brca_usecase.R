library(dplyr)
library(tidyverse)
library(ggplot2)

setwd("/home/lmgonzalezb/Documents/bc-meta/cohort_statistics/")
source("../bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")

subclones <- read.table("subclone_level_annotated.tsv", header = TRUE)


theme_barplot <- theme(
  panel.background = element_blank(),
  plot.background = element_blank(), 
  panel.grid.major = element_line(color = "grey90"),
  panel.grid.minor = element_blank(),
  plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
  axis.title.x = element_text(size = 14, margin = margin(t = 6), face = "bold"),
  axis.title.y = element_text(size = 14, margin = margin(r = 6), face = "bold"),
  axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
  axis.text.y = element_text(size = 12, color = "black"),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12)
)

#### Focus on BRCA ###
subclones_brca <- subclones %>%
  mutate(Study_Sample = paste0(Study, "_", Sample)) %>%
  filter(Refined.Tumor.Type == "BRCA") %>%
  select(-c(2:6)) %>%
  mutate(
    scTherapy.Cluster = factor(scTherapy.Cluster, levels = as.character(1:10)),
    Treated = ifelse(Treated, "Treated", "Untreated"),
    Sex = ifelse(Sex == "f", "Female", "Male"),
    Sample.Type = factor(ifelse(Sample.Type == "m", "Metastasis", "Primary"),
                         levels = c("Primary", "Metastasis")),
    Age.Group = factor(cut(Age, breaks = c(0, 45, 55, Inf), right = FALSE,
                    labels = c("<45", "45-55", "+55")),
                    levels =  c("<45", "45-55", "+55"))
  ) %>%
  mutate(across(everything(), ~ if_else(. == "", "Unknown", .))) %>%
  distinct()

subclones_brca[is.na(subclones_brca)] <- 
  
subclones_brca_long <- subclones_brca %>%
  select(Sample.Type, TME.Archetype, Tumor.Subtype, Age.Group, scTherapy.Cluster) %>%
  pivot_longer(cols = c(Sample.Type, TME.Archetype, Tumor.Subtype, Age.Group),
               names_to = "Variable",
               values_to = "Category") %>%
  mutate(Category = factor(Category, levels = c("<45", "45-55", "+55", 
                                                "Primary", "Metastasis",
                                                sort(TME.Archetype),
                                                )))


barplot <- ggplot(subclones_brca_long, aes(x = Category,
                                      fill = scTherapy.Cluster)) +
  geom_bar(position = "fill", width = 0.7) +
  facet_grid(~ Variable, scales = "free_x", space = "free_x") +
  labs(
    y = "Fraction of subclones",
    fill = "scTherapy cluster"
  ) +
  scale_fill_manual(values = sctherapy_colors) +
  theme_barplot +
  theme(panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1))

subclones_pct <- subclones_brca_long %>%
  group_by(Variable, Category, scTherapy.Cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Variable, Category) %>%
  mutate(pct = count / sum(count)) %>%
  ungroup()

# Graficar barras agrupadas con proporciones
ggplot(subclones_pct, aes(x = Category, y = pct, fill = scTherapy.Cluster)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(~ Variable, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Percentage", fill = "scTherapy cluster") +
  scale_fill_manual(values = sctherapy_colors) +
  theme_barplot +
  theme(
    panel.spacing = unit(1.5, "lines"),
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1)
  )
ggsave("./figures/treatede_cluster_tumortype.png", plot = barplot, width = 16, height = 10)


barplot <- ggplot(subclones_brca_long, aes(x = Category,
                                           fill = scTherapy.Cluster)) +
  geom_bar(position = "stack", width = 0.7) +
  facet_grid(~ Variable, scales = "free_x", space = "free_x") +
  labs(
    y = "Fraction of subclones",
    fill = "scTherapy cluster"
  ) +
  scale_fill_manual(values = sctherapy_colors) +
  theme_barplot +
  theme(panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1))


samples_brca <- subclones_brca %>%
  select(c(3:27)) %>%
  distinct()

barplot <- ggplot(samples_brca, aes(x = Age.Group,
                                           fill = TME.Archetype)) +
  geom_bar(position = "fill", width = 0.7) +
  labs(
    y = "Fraction of subclones",
    fill = "TME archetype"
  ) +
  scale_fill_manual(values = tme_colors) +
  theme_barplot +
  theme(panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1))


## MoAs of drugs predicted for subclones within each age
subclones_brca_MoAs <- subclones %>%
  mutate(Study_Sample = paste0(Study, "_", Sample)) %>%
  filter(Refined.Tumor.Type == "BRCA") %>%
  mutate(
    Age.Group = factor(cut(Age, breaks = c(0, 45, 55, Inf), right = FALSE,
                           labels = c("<45", "45-55", "+55")),
                       levels =  c("<45", "45-55", "+55"))
  ) %>%
  select(1,2,4, Age.Group) %>%
  distinct()


barplot <- ggplot(subclones_brca_MoAs, aes(x = Age.Group,
                                       fill = Drug.Mechanism.of.Action)) +
  geom_bar(position = "fill", width = 0.7) +
  labs(
    y = "Fraction of subclones",
    fill = "Age group"
  ) +
  scale_fill_manual(values = MoAs_colors) +
  theme_barplot +
  guides(fill = guide_legend(ncol = 1))


