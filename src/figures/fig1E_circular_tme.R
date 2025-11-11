library(dplyr)
library(tidyverse)
library(ggplot2)

setwd("/home/lmgonzalezb/Documents/bc-meta/cohort_statistics/")
source("../bc-meta_repo/bc-meta/src/figures/TCCA_palette.R")

theme_barplot <- theme(
  plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
  axis.title.x = element_text(size = 14, margin = margin(t = 6), face = "bold"),
  axis.title.y = element_text(size = 14, margin = margin(r = 6), face = "bold"),
  axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
  axis.text.y = element_text(size = 12, color = "black"),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12)
)

metadata <- read.table("tcca_metadata.tsv", sep = "\t", header = TRUE)

## Create a concentric plot with sample fraction for general and specific TME archetypes
samples_with_tme <- metadata %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  group_by(study_sample) %>%
  summarise(has_true = any(malignancy == "True"),
            has_false = any(malignancy == "False")) %>%
  filter(has_true & has_false) %>%
  pull(study_sample)

metadata_with_tme <- metadata %>%
  mutate(study_sample = paste0(study, "_", sample)) %>%
  filter(study_sample %in% samples_with_tme & !(refined_tumor_type %in% c("ALL", "CLL", "LAML", "MM")))

tme_subtypes <- metadata_with_tme %>%
  mutate(
    tme_split = str_split_fixed(tme_archetype, "_", n = 3),
    `tme_general` = if_else(
      tme_split[,3] == "", tme_archetype,
      paste(tme_split[,1], tme_split[,2], sep = "_")
    ),
    `tme_specific` = if_else(
      tme_split[,3] == "", tme_archetype,
      str_remove(tme_archetype, paste0("^", tme_split[,1], "_", tme_split[,2], "_"))
    )
  ) %>%
  select(-tme_split) %>%
  select(study, sample, refined_tumor_type, tme_general, tme_specific) %>%
  distinct()


df <- metadata_with_tme %>%
  count(tme_general, tme_specific, name = "count") %>%
  group_by(tme_general) %>%
  mutate(
    general_total = sum(count),
    specific_frac = count / general_total
  ) %>%
  ungroup() %>%
  mutate(
    overall_total = sum(count),
    general_frac = general_total / overall_total
  )

df <- df %>%
  group_by(tme_general, tme_specific) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  ungroup()

total <- sum(df$count)

# Fractions for plotting
df <- df %>%
  arrange(tme_general, tme_specific) %>%
  mutate(
    specific_frac = count / total,
    specific_ymax = cumsum(specific_frac),
    specific_ymin = lag(specific_ymax, default = 0)
  )

# Now compute general arcs (outer ring)
df_general <- df %>%
  group_by(tme_general) %>%
  summarise(
    general_frac = sum(specific_frac),
    general_ymax = max(specific_ymax),
    general_ymin = min(specific_ymin)
  )

tme_circular <- ggplot() +
  # Center sectors for tme_specific
  geom_rect(data = df,
            aes(
              ymin = specific_ymin, ymax = specific_ymax,
              xmin = 0.3, xmax = 1,  # inner circle width
              fill = tme_specific
            ),
            color = "white") +
  
  # Outer donut for tme_general
  geom_rect(data = df_general,
            aes(
              ymin = general_ymin, ymax = general_ymax,
              xmin = 1, xmax = 1.4,  # outer ring
              fill = tme_general
            ),
            color = "white") +
  
  coord_polar(theta = "y") +
  xlim(0, 1.5) +
  theme_void() +
  theme(legend.position = "right") +
  ggtitle("TME Specific Sector Plot with Outer TME General Ring")

ggsave("figures/tme_circular_plot.pdf", plot = tme_circular, width = 10, height = 10)




### Compute TME distributions depending on sex and age groups
clinical <- metadata_with_tme  %>%
  select(3, 9:29, 33) %>%
  distinct() %>%
  mutate(
    age_group = cut(
      age,
      breaks = c(seq(0, 100, by = 10), Inf),  # Final bin: 100+
      include.lowest = TRUE,
      right = FALSE,
      labels = c(
        paste(seq(0, 90, by = 10), seq(10, 100, by = 10), sep = "-"),
        "100+"
      )
    ),
    sex = ifelse(sex == "", "Unknown", sex),
    age_group = ifelse(is.na(as.character(age_group)), "Unknown", as.character(age_group)),
    tme_archetype = factor(tme_archetype, 
                           levels = c(setdiff(sort(unique(tme_archetype)), 
                                              "none"), "none"))
  )

# Compute number of samples within each sex and age group 
counts <- clinical %>%
  group_by(age_group, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(sex_label = paste0(sex, "\n(n=", n, ")"))

# Step 2: Merge back to clinical data to use custom labels
clinical_labeled <- clinical %>%
  left_join(counts, by = c("age_group", "sex")) %>%
  mutate(sex_label = factor(sex_label, levels = unique(sex_label)))

# Step 3: Plot using sex_label on the x-axis
barplot <- ggplot(clinical_labeled, aes(x = sex_label, fill = tme_archetype)) +
  geom_bar(position = "fill") +
  facet_grid(. ~ age_group, scales = "free_x", space = "free_x") +
  labs(
    title = "TME Archetype Counts by Age Group and Sex",
    x = "Age group",
    y = "Proportion of samples",
    fill = "TME archetype"
  ) +
  scale_fill_manual(values = tme_colors) +
  scale_x_discrete(drop = TRUE) +
  theme_barplot +
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(size = 14)
  ) +
  guides(fill = guide_legend(ncol = 3))

barplot
ggsave(
  "./figures/tme/tme_across_sex_age.pdf",
  plot = barplot,
  width = 15,
  height = 8,
  dpi = 300
)
