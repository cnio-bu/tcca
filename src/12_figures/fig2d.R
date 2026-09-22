library(GenVisR)
library(tidyverse)
library(ggplot2)
library(GenomicRanges)
library(patchwork)
library(AnnotationHub)
library(dplyr)
library(ComplexHeatmap)

# Load CNV data
setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/cna_metadata")
cs <- read.table("cnv_segments_clones_lvl2_cytobands.tsv", header = T, sep = "\t")

# This is a table of averages of amplification and deletions of each region
amp <- apply(cs, 1, function(row) mean(row[row > 2], na.rm = TRUE))
del <- apply(cs, 1, function(row) mean(row[row < 2], na.rm = TRUE))
amp.del <- cbind(amp, del)
amp.del <- replace(amp.del, is.nan(amp.del), 0)

# Make a granges of the cytobands of the human genome
proxy <- "mgonzalezb@cnio.es"
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

# Join cytobands to CNV data and set long format
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


# Load clinical and TCs metadata
metadata <- read.table("/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/tcca_metadata_h5ad.tsv", header = T, sep = "\t")
metadata <- metadata %>%
    mutate(study_sample = paste(study, sample, sep = "__")) %>% 
    select(scevan_subclone, study, study_sample, tumor_type, sample_type, treated, therapeutic_cluster) %>%
    filter(!is.na(therapeutic_cluster)) %>%
    dplyr::rename(sample = scevan_subclone) %>%
    distinct()


# Add samples IDs and other metadata to cs.long
cs.long <- cs.long %>%
    mutate(sample = gsub("subclone", "", str_replace_all(sample, "__", "\\.")))
cs.long.merged <- left_join(cs.long, metadata, by = "sample")


# Plot heatmap of CNV data per TC
ht_data <- cs.long.merged %>%
  filter(study != "cell_lines_gabriella_kinker")

cnv_freq_list <- list()

for (tc in 1:10) {
    tc_samples <- ht_data %>% filter(therapeutic_cluster == tc)
    tc_freq <- cnFreq(tc_samples,
                      CN_low_cutoff = 1.5,
                      CN_high_cutoff = 2.5,
                      genome = "hg38",
                      out = "data")
    tc_freq$data$TC <- paste0("TC", tc)
    cnv_freq_list[[paste0("TC", tc)]] <- tc_freq$data
}

# Combine frequencies for all
cnv_freq_all <- bind_rows(cnv_freq_list)


# Add region name
region_name <- ht_data %>%
    select(region_name, chromosome, start, end) %>%
    distinct()

cnv_freq_all <- cnv_freq_all %>%
    left_join(region_name, by = c("chromosome", "start", "end")) %>%
    mutate(chr_numeric = as.integer(gsub("chr", "", chromosome)))

region_order <- cnv_freq_all %>%
    distinct(region_name, chromosome, start, end, chr_numeric) %>%
    arrange(chr_numeric, start) %>%
    pull(region_name)

# Gain matrix (TCs x regions)
gain_matrix <- cnv_freq_all %>%
    dplyr::select(TC, region_name, gainProportion) %>%
    pivot_wider(
        names_from = region_name,
        values_from = gainProportion,
        values_fill = 0
    ) %>%
    column_to_rownames("TC") %>%
    as.matrix()

gain_matrix <- gain_matrix[paste0("TC", 1:10), region_order]

# Loss matrix
loss_matrix <- cnv_freq_all %>%
    dplyr::select(TC, region_name, lossProportion) %>%
    pivot_wider(
        names_from = region_name,
        values_from = lossProportion,
        values_fill = 0
    ) %>%
    column_to_rownames("TC") %>%
    as.matrix()

loss_matrix <- loss_matrix[paste0("TC", 1:10), region_order]

# Combined matrix: gains are positive, losses are negative
combined_matrix <- gain_matrix - loss_matrix

# Chromosome top annotation
region_chr <- cnv_freq_all %>%
    distinct(region_name, chr_numeric) %>%
    arrange(match(region_name, region_order)) %>%
    pull(chr_numeric)

block_colors <- ifelse(unique(region_chr) %% 2 == 0, "#bdbdbd", "#e2e2e2")
block_labels <- 1:22

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    foo = anno_block(
        gp = gpar(fill = block_colors),
        labels = block_labels,
        labels_gp = gpar(fontsize = 8)
    )
)

# Cytobands of interest to highlight
regions_interest <- c(
    "1q21.3", "12p13.2", "12q13.13", "17q21.2", "10p12.31", 
    "17q21.1", "13q13.3", "13q14.11", "7p15.1", "7p21.1", "3q25.1", "8q24.3",
    "11q13.2", "19q13.2", "2q31.3", "3p11.1", "3q11.1", "5q14.1",
    "10p12.1", "13q12.11", "13q32.1", "21q22.3", "7p22.1", "7q31.1",
    "1q23.3", "6q23.2", "10p11.1", "10q11.1", "11q12.1", "17q11.1",
    "22q13.31"
)
regions_interest <- paste0("chr", regions_interest)
mark_idx <- which(region_order %in% regions_interest)
mark_labels <- region_order[mark_idx]

# Bottom annotation: mark regions of interest
bottom_annotation <- ComplexHeatmap::HeatmapAnnotation(
    mark = anno_mark(
        at = mark_idx,
        labels = mark_labels,
        labels_gp = gpar(fontsize = 7),
        side = "bottom",
        link_height = unit(3, "mm")
    ),
    which = "column"
)


# Heatmap color
cnv_colors <- circlize::colorRamp2(
  breaks = seq(-0.5, 0.5, length.out = 11),
  colors = rev(RColorBrewer::brewer.pal(11, "RdBu"))
)


# Plot heatmap
hm_cnv <- Heatmap(
    matrix = combined_matrix,
    name = "CNV\nproportion",
    col = cnv_colors,
    top_annotation = top_annotation,
    bottom_annotation = bottom_annotation,
    
    # Rows (TCs)
    cluster_rows = FALSE,
    row_order = paste0("TC", 1:10),
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    
    # Columns (regions)
    cluster_columns = FALSE,
    column_order = region_order,
    show_column_names = FALSE,
    column_split = region_chr,
    column_gap = unit(0.5, "mm"),
    column_title = NULL,
    
    # Legend
    heatmap_legend_param = list(
        title = "Net CNA\nproportion\n(gain - loss)",
        title_gp = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 8),
        direction = "vertical"
    ),
    
    width = unit(20, "cm"),
    height = unit(6, "cm")
)

pdf(
    file = "/storage/scratch01/shared/projects/bc-meta/single_cell/seurat/tcca/heatmap_cnv_TC.pdf",
    #res = 300,
    width = 14,
    height = 5,
    #units = "in"
)

draw(hm_cnv)

dev.off()