library(BPCells)
library(ComplexHeatmap)
library(Seurat)
library(tidyverse)


setwd("/local/bc_meta/")
source(file = "/local/bc_meta/TCCA_palette.R")

## Load sketched matrix and metadata

sketched_mat <- readRDS("v5/full_genes_copynumber_sketch_5k_scaled_normcounts_mtx_lvl2.rds")

metadata <- read.table("v5/sketch_mat_metadata_5k_lvl2.tsv", header = T)
metadata <- metadata %>%
  mutate(
    barcode_study_sample = paste(original_barcode, study, sample, sep = "__")
  ) %>%
  filter(!(is.na(malignancy)))

sketched_mat <- sketched_mat[, rownames(metadata)]

## Get gene info df
original_rownames <- rownames(sketched_mat)

split_names <- strsplit(original_rownames, "-")
split_names_df <- as.data.frame(do.call(rbind, split_names), stringsAsFactors = FALSE)

## LRG genes contain a "-" so they need special preprocessing
lrg <- filter(split_names_df, V4 == "LRG") %>%
  unite(newV4, V4, V5, sep = "-")
colnames(lrg) <- c("chromosome", "start", "end", "name", "symbol")

not_lrg <- filter(split_names_df, V4 != "LRG")
colnames(not_lrg) <- c("chromosome", "start", "end", "name", "symbol")

split_names_df <- bind_rows(lrg, not_lrg[, 1:5]) %>%
  mutate(
    chromosome_numeric = as.integer(sub("X", "", chromosome)),
    start = as.numeric(start),
    original_name = paste(chromosome, start, end, name, symbol, sep = "-")
  )

## Order genes by chromosomic position and reorder matrix
split_names_df <- split_names_df %>%
  arrange(chromosome_numeric, start)

reorder_idx <- match(split_names_df$original_name, rownames(sketched_mat))
sketched_mat <- sketched_mat[reorder_idx, ]

## Assign gene symbols to rownames
rownames(sketched_mat) <- split_names_df$symbol
sketched_mat <- t(sketched_mat)

## Bin matrix to reduce dimensionality in genes

grouped_list <- list()
grouped_names <- c()

stopifnot(all(colnames(sketched_mat) == split_names_df$symbol))

chromosomes <- unique(split_names_df$chromosome_numeric)

for (chr in chromosomes) {
  chr_idx <- which(split_names_df$chromosome_numeric == chr)
  chr_mat <- sketched_mat[, chr_idx, drop = FALSE]
  chr_names <- colnames(chr_mat)
  
  for (i in seq(1, ncol(chr_mat), by = 3)) {
    end_idx <- min(i + 2, ncol(chr_mat))
    cols <- chr_mat[, i:end_idx, drop = FALSE]
    
    if (ncol(cols) == 1) {
      group_mean <- cols[, 1]
    } else {
      group_mean <- rowMeans(cols)
    }
    
    grouped_list[[length(grouped_list) + 1]] <- group_mean
    
    name_start <- chr_names[i]
    name_end <- chr_names[end_idx]
    grouped_names <- c(grouped_names, paste0(name_start, "_to_", name_end))
  }
}

grouped_mat <- do.call(cbind, grouped_list)
colnames(grouped_mat) <- grouped_names

## Remove some malignant cells with few CNV for plotting
row_sd <- apply(grouped_mat, 1, sd)

malignant <- metadata$malignancy == TRUE
rows_to_remove <- which(row_sd < 0.075 & malignant)

grouped_mat <- grouped_mat[-rows_to_remove, ]
metadata <- metadata[-rows_to_remove, ]

## Set heatmap annotation
## Clones annotation
clones_annot <- metadata %>% dplyr::select(malignancy) %>%
  mutate(malignancy = ifelse(malignancy == TRUE, "Malignant", "Healthy"))

pals <- list(
  malignancy = c(
    Malignant = "#db4646",
    Healthy = "#0090ab"
  )
)

left_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df =  clones_annot,
  which = "row",
  col = pals
)

## Gene groups annotations
grouped_symbols <- sapply(strsplit(colnames(grouped_mat), "_to_"), `[`, 1)

genes_annot <- split_names_df %>%
  filter(symbol %in% grouped_symbols) %>%
  distinct(symbol, chromosome_numeric) %>%
  column_to_rownames("symbol") %>%
  dplyr::select(chromosome_numeric)

genes_annot <- genes_annot[grouped_symbols, , drop = FALSE]
colnames(genes_annot) <- "Chromosome"
rownames(genes_annot) <- colnames(grouped_mat)

## Assign colors: dark gray for odd, light gray for even
block_labels <- 1:22
block_colors <- ifelse(block_labels %% 2 == 0, "#bdbdbd", "#e2e2e2")

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  foo = anno_block(gp = gpar(fill = block_colors), labels = block_labels)
)

hm <- ComplexHeatmap::Heatmap(
  mat = grouped_mat,
  name = "Segm.mean",
  left_annotation = left_annotation,
  top_annotation = top_annotation,
  row_split = clones_annot$malignancy,
  
  # Rows
  cluster_rows = FALSE,
  row_order = rownames(grouped_mat),
  row_title = "Cells",
  show_row_names = FALSE,
  
  # Cols
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_order = colnames(grouped_mat),
  column_title = "Chromosomes",
  column_split = genes_annot$Chromosome,
  
  heatmap_width = unit(6, "in"),
  heatmap_height = unit(10, "in")
)

png(
  file = "heatmap_tumor_healthy_cells.png",
  res = 300,
  width = 10,
  height = 20,
  units = "in"
)

draw(hm)
dev.off()
