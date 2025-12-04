library(GenVisR)
library(tidyverse)
library(ggplot2)
library(GenomicRanges)
library(patchwork)
library(AnnotationHub)

## Load CNV data
setwd("/local/bc_meta/scevan/cytobands_CNV/")
cs <- read.table("/local/bc_meta/scevan/cnv_segments_clones_lvl2_cytobands.tsv", header = T, sep = "\t")

## This is a table of averages of amplification and deletions of each region
amp <- apply(cs, 1, function(row) mean(row[row > 2], na.rm = TRUE))
del <- apply(cs, 1, function(row) mean(row[row < 2], na.rm = TRUE))
amp.del <- cbind(amp, del)
amp.del <- replace(amp.del, is.nan(amp.del), 0)

## Make a granges of the cytobands of the human genome
proxy <- "lserranor@cnio.es"
proxy <- httr::use_proxy(Sys.getenv('http_proxy'))
httr::set_config(proxy)
AnnotationHub::setAnnotationHubOption("PROXY", proxy)
AnnotationHub::getAnnotationHubOption("LOCAL")

hub <- AnnotationHub()

hub_hg38 <- subset(hub, 
                   (hub$species == "Homo sapiens") & (hub$genome == "hg38"))

cytobands  <- hub_hg38[[797]]
cytobands$custom_name <- paste0(seqnames(cytobands), cytobands$name)

dfA <- data.frame(names = cytobands$custom_name,
                  chrom = as.data.frame(seqnames(cytobands)))
dfB <- as.data.frame(ranges(cytobands))

cn.data <- cbind(dfA, dfB)
cn.data <- cn.data[,-5]

## Join cytobands to CNV data and set long format
cs <- rownames_to_column(as.data.frame(cs), var = "names")
cs <- merge(cn.data, cs, by = "names")

cs.long <- cs %>%
  pivot_longer(-c(names, value, start, end), names_to = "sampleID", values_to = "mean")

cs.long <- data.frame(region_name = cs.long$names,
                      sample = cs.long$sampleID,
                      chromosome = cs.long$value,
                      start = cs.long$start,
                      end = cs.long$end,
                      probes = 1,
                      segmean = cs.long$mean)

cs.long <- replace(cs.long, is.na(cs.long), 2) %>%
  mutate(sample = str_replace_all(sample, "\\.", "-")) # This reverts scevan name changing

cs.long$sample <- gsub("GSM5645908_Breast_1_biol-rep", "GSM5645908_Breast_1_biol.rep", cs.long$sample) # This three fix three exceptions
cs.long$sample <- gsub("SyS11-met", "SyS11.met", cs.long$sample)
cs.long$sample <- gsub("Travaglini_Krasnow_2020_distal-1b", "Travaglini_Krasnow_2020_distal 1b", cs.long$sample)


## Calculate ITH 
samples <- strsplit(unique(cs.long$sample), "__")
samples <- as.data.frame(do.call(rbind, samples))
rownames(samples) <- unique(cs.long$sample)
samples <- rownames_to_column(samples, var = "sample")
colnames(samples) <- c("sample","study", "sample_original", "clone")
samples <- mutate(samples,
                  study__sample = paste(study, sample_original, sep = "__"))

## Load clinical and TCs metadata
clinical <- read.table("/local/bc_meta/clinical_metadata_v4_clean.tsv", header = T, sep = "\t")
clinical <- clinical %>%
  mutate(true_patient = paste0(study, "_", patient),
         study__sample = paste(study, sample, sep = "__")) %>%
  filter(study__sample %in% samples$study__sample) # Keep samples for which we have CNV data

metacom <- readRDS("/local/bc_meta/scevan/threapeutic_clusters.rds")
metacom <- metacom %>%
  mutate(sample = str_replace_all(subclone, "\\.", "__") %>% str_replace("(\\d+)$", "subclone\\1"))
metacom$sample <- gsub("GSM5645908_Breast_1_biol__rep", "GSM5645908_Breast_1_biol.rep", metacom$sample) # Fix manual exceptions
metacom$sample <- gsub("SyS11__met", "SyS11.met", metacom$sample)

## Add samples IDs and other metadata to cs.long
cs.long.merged <- left_join(cs.long, samples, by = "sample")
cs.long.merged <- left_join(cs.long.merged, metacom, by = "sample")
cs.long.merged <- left_join(cs.long.merged, 
                            dplyr::select(clinical, study__sample, tumor_type, stage, sample_type, treated, refined_tumor_site), 
                            by = "study__sample")


sigs_to_keep <- c("chr3q25.1", "chr3q25.2", "chr3q26.32", "chr3q26.33", 
                  "chr3q27.1", "chr3q27.2", "chr3q27.3", "chr3q28", "chr3q29")

data0 <- cs.long.merged %>%
  filter(cluster == 4,
         region_name %in% sigs_to_keep
  )

data0  %>%
  distinct(sample_original, tumor_type) %>%  
  count(tumor_type, name = "n_unique_samples")

all_combo_counts <- list()

for (cancer in unique(data0$tumor_type)) {
  all_samples <- data0 %>%
    filter(
      tumor_type == cancer
    ) %>%
    distinct(sample)
  
  data <- data0 %>%
    filter(
      tumor_type == cancer
    ) %>%
  mutate(is.amplified = ifelse(segmean >= 2.5, TRUE, FALSE)) %>%
  dplyr::select(region_name, sample, segmean, is.amplified)
  
  total_n <- nrow(all_samples)  
  
  # Get amplified combos by sample
  amplified_combos <- data %>%
    filter(is.amplified) %>%
    group_by(sample) %>%
    summarise(amplified_combo = paste(sort(unique(region_name)), collapse = ","), .groups = "drop")
  
  # Assign "" to samples without any amplification
  combo_per_sample <- all_samples %>%
    left_join(amplified_combos, by = "sample") %>%
    mutate(amplified_combo = ifelse(is.na(amplified_combo), "", amplified_combo))
  
  n_amplified <- combo_per_sample %>%
    filter(amplified_combo != "") %>%
    nrow()
  
  # Count combos and %
  combo_counts <- combo_per_sample %>%
    group_by(amplified_combo) %>%
    summarise(n_samples = n(), .groups = "drop") %>%
    mutate(
      percentage_total = 100 * n_samples / total_n,
      percentage_amplified = ifelse(amplified_combo != "", 100 * n_samples / n_amplified, NA)
    ) %>%
    arrange(desc(percentage_total))
  
  all_combo_counts[[cancer]] <- combo_counts
}

## Without filtering by cancer type
all_samples <- data0 %>%
  distinct(sample)

data <- data0 %>%
  mutate(is.amplified = ifelse(segmean >= 2.5, TRUE, FALSE)) %>%
  dplyr::select(region_name, sample, segmean, is.amplified)

total_n <- nrow(all_samples)  

# Get amplified combos by sample
amplified_combos <- data %>%
  filter(is.amplified) %>%
  group_by(sample) %>%
  summarise(amplified_combo = paste(sort(unique(region_name)), collapse = ","), .groups = "drop")

# Assign "" to samples without any amplification
combo_per_sample <- all_samples %>%
  left_join(amplified_combos, by = "sample") %>%
  mutate(amplified_combo = ifelse(is.na(amplified_combo), "", amplified_combo))

n_amplified <- combo_per_sample %>%
  filter(amplified_combo != "") %>%
  nrow()

# Count combos and %
combo_counts <- combo_per_sample %>%
  group_by(amplified_combo) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  mutate(
    percentage_total = 100 * n_samples / total_n,
    percentage_amplified = ifelse(amplified_combo != "", 100 * n_samples / n_amplified, NA)
  ) %>%
  arrange(desc(percentage_total))

all_combo_counts[["All"]] <- combo_counts


library(ComplexUpset)
source("/local/bc_meta/TCCA_palette.R")

# Extraer todos los combos en una tabla long
binary_df <- combo_per_sample %>%
  mutate(region_list = strsplit(amplified_combo, ",")) %>%
  unnest(region_list, keep_empty = TRUE) %>%
  mutate(value = TRUE) %>%
  pivot_wider(
    id_cols = sample,
    names_from = region_list,
    values_from = value,
    values_fill = FALSE
  )

tumor_type_data <- data0 %>% dplyr::select(sample, tumor_type) %>% unique()
binary_df <- left_join(binary_df, 
                       tumor_type_data,
                       by = "sample")

cnv_order <- c("chr3q25.1", "chr3q25.2", "chr3q26.32", "chr3q26.33",
               "chr3q27.1", "chr3q27.2", "chr3q27.3", "chr3q28", "chr3q29")


ComplexUpset::upset(
  binary_df,
  intersect = cnv_order,
  mode = "exclusive_intersection",
  name = "Amplified regions",
  width_ratio = 0.15,
  sort_sets=FALSE,
  base_annotations = list(
    'Intersection size' = ComplexUpset::intersection_size(
      mode = "exclusive_intersection",
      counts = FALSE,
      mapping = aes(fill = tumor_type)
    ) +
      scale_fill_manual(values = tumor_type_colors)
    ),
  guides = "over",
  set_sizes=(
    upset_set_size()
    + theme(axis.text.x=element_text(angle=90))
  )
)

ggsave("/local/bc_meta/gdsc/upset_plot_cnvs_tc4.png", width = 8, height = 7, dpi = 300)
