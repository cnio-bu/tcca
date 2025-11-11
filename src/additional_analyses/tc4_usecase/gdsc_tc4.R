library(BPCells)
library(Seurat)
library(tidyverse)
library(GenVisR)
library(ggplot2)
library(GenomicRanges)
library(patchwork)
library(AnnotationHub)
library(GSEABase)
library(UCell)
library(dplyr)
library(purrr)
library(broom)

setwd("/local/bc_meta/")
source("TCCA_palette.R")


library(readxl)
library(stringr)

# Leer el archivo Excel
df <- read_excel("/home/lserranor/Downloads/top_drugs_cluster.xlsx")

# Función de formato personalizada
# Función de formato personalizada que maneja NA
formatear_nombre <- function(x) {
  sapply(x, function(drug) {
    if (is.na(drug)) {
      return(NA)
    } else if (str_detect(drug, "\\d")) {
      return(toupper(drug))  # Si contiene número, todo en mayúsculas
    } else {
      return(str_to_title(tolower(drug)))  # Si no, primera letra mayúscula
    }
  }, USE.NAMES = FALSE)
}
# Crear y formatear los vectores
cluster_1  <- formatear_nombre(df$Cluster_1)
cluster_2  <- formatear_nombre(df$Cluster_2)
cluster_3  <- formatear_nombre(df$Cluster_3)
cluster_4  <- formatear_nombre(df$Cluster_4)
cluster_5  <- formatear_nombre(df$Cluster_5)
cluster_6  <- formatear_nombre(df$Cluster_6)
cluster_7  <- formatear_nombre(df$Cluster_7)
cluster_8  <- formatear_nombre(df$Cluster_8)
cluster_9  <- formatear_nombre(df$Cluster_9)
cluster_10 <- formatear_nombre(df$Cluster_10)

all_drugs <- c(
  cluster_1, cluster_2, cluster_3, cluster_4, cluster_5,
  cluster_6, cluster_7, cluster_8, cluster_9, cluster_10
)

## FUNCTIONS -------------------------

## Extract desired strudy from global seurat object 
metadata <- read.table("tcca_metadata.tsv", header = T, sep = "\t") %>%
  column_to_rownames("cell")

extract_seu <- function(path, sam = NULL) {
  if (is.null(sam)) {
    mat <- open_matrix_dir(path)
    seu <- CreateSeuratObject(mat,
                              meta.data = metadata)
  }
  
  else {
    mat <- open_matrix_dir(path)
    seu <- CreateSeuratObject(mat,
                              meta.data = metadata)
    seu <- subset(seu, sample %in% sam)
  }
  
  return(seu)
}

## Read GMT as list
read_gmt_list <- function(gmt) {
  gmt_lines <- readLines(gmt)
  
  signatures <- lapply(gmt_lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    genes <- parts[-c(1, 2)]  # Remove name and description
    return(genes)
  })
  
  # Name signatures
  names(signatures) <- sapply(gmt_lines, function(line) strsplit(line, "\t")[[1]][1])
  
  return(signatures)
}

## Plot correlation from final data table
plot_correlation <- function(data, drug, tumor_type = "all", cluster) {
  
  if (unique(tumor_type == "all")) {
    df_plot <- data %>%
      filter(DRUG_NAME == drug) %>%
      dplyr::select(all_of(cluster), AUC, refined_tumor_type) %>%
      na.omit()
  } else {
    df_plot <- data %>%
      filter(DRUG_NAME == drug,
             refined_tumor_type %in% tumor_type
      ) %>%
      dplyr::select(all_of(cluster), AUC, refined_tumor_type) %>%
      na.omit()
  }
  
  if (nrow(df_plot) < 3) {
    stop("Not enough samples (3) to calculate correlation")
  }
  
  # Correlation test
  test <- cor.test(df_plot[[cluster]], df_plot$AUC, method = "pearson")
  corr_val <- round(test$estimate, 3)
  p_val <- test$p.value
  p_val_formatted <- ifelse(p_val < 0.001, "< 0.001", paste0("= ", round(p_val, 4)))
  
  # Plot
  if (unique(tumor_type == "all")) {
    ggplot(df_plot, aes(x = .data[[cluster]], y = AUC, color = refined_tumor_type)) +
      geom_point(size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      labs(
        title = paste0("Pearson correlation: ", cluster, " and AUC"),
        subtitle = paste0("Drug: ", drug, " | All tumor types", "\n", "Pearson r = ", corr_val, ", p ", p_val_formatted),
        x = cluster,
        y = "AUC",
        color = "Tumor_type"
      ) +
      scale_color_manual(values = tumor_type_colors) +
      theme_minimal()
  } else {
    tumor_text <- paste(tumor_type, collapse = ", ")
    ggplot(df_plot, aes(x = .data[[cluster]], y = AUC, color = refined_tumor_type)) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      labs(
        title = paste0("Pearson correlation: ", cluster, " and AUC"),
        subtitle = paste0("Drug: ", drug, " | ", tumor_text, "\n", "Pearson r = ", corr_val, ", p ", p_val_formatted),
        x = cluster,
        y = "AUC",
        color = "Tumor_type"
      ) +
      scale_color_manual(values = tumor_type_colors) +
      theme_minimal()
  }
}

## -----------------------

## Get Seurat objects and GDSC2 data
gdsc <- read.table("gdsc/GDSC2_fitted_dose_response_27Oct23.tsv", sep = "\t", header = T) %>%
  mutate(CELL_LINE_NAME = sub("-", "", CELL_LINE_NAME)) %>%
  dplyr::select(CELL_LINE_NAME, DRUG_NAME, PATHWAY_NAME, AUC) %>%
  group_by(CELL_LINE_NAME, DRUG_NAME, PATHWAY_NAME) %>%
  summarise(AUC = mean(AUC), .groups = "drop") ## These steps averages duplicated AUC (same drug, different ID)

seu <- extract_seu("v5/lvl2/cell_lines_gabriella_kinker_v5/")

seu@meta.data <- seu@meta.data %>%
  mutate(CELL_LINE_NAME = sub("_.*", "", sample))

common_cell_lines <- intersect(seu$CELL_LINE_NAME, gdsc$CELL_LINE_NAME)

seu <- subset(seu, subset = CELL_LINE_NAME %in% common_cell_lines)
gdsc <- subset(gdsc, subset = CELL_LINE_NAME %in% common_cell_lines)


markers <- read_gmt_list("gdsc/tc4_cnvs.gmt")

## Adapt signatures to keep only genes expressed in the dataset

expr_matrix <- GetAssayData(seu, slot = "counts")  # o "data" para valores normalizados
expr_genes <- rownames(expr_matrix)[Matrix::rowSums(expr_matrix > 0) > 0]

markers <- lapply(markers, function(genes_vector) {
  intersect(genes_vector, expr_genes)
})


seu <- AddModuleScore_UCell(seu, features=markers, name = NULL)

## Get dataframe of averages

tc_sigs <- seu@meta.data %>%
  dplyr::select(CELL_LINE_NAME, refined_tumor_type, matches("^[0-9]{1,2}q")) %>%
  group_by(CELL_LINE_NAME, refined_tumor_type) %>%
  summarise(across(matches("^[0-9]{1,2}q"), mean, na.rm = TRUE), .groups = "drop")

data <- left_join(tc_sigs, gdsc)


tc4_drugs <- intersect(cluster_4, data$DRUG_NAME)

## Calculate correlations by combination DRUG_NAME + CELL_LINE_NAME
cluster_cols <- grep("^[0-9]{1,2}q", colnames(data), value = TRUE)

results <- data %>%
  group_by(DRUG_NAME, refined_tumor_type) %>%
  group_map(~ {
    drug <- unique(.x$DRUG_NAME)
    cell <- unique(.x$refined_tumor_type)
    moa <- unique(.x$PATHWAY_NAME)
    
    ## Check for enough replicates
    if (nrow(.x) < 3) return(NULL)
    
    map_dfr(cluster_cols, function(cluster) {
      # Error control if NA o variance = 0
      tryCatch({
        test <- cor.test(.x[[cluster]], .x$AUC, method = "pearson")
        tibble(
          DRUG_NAME = drug,
          PATHWAY_NAME = moa,
          refined_tumor_type = cell,
          Cluster = cluster,
          Correlation = test$estimate,
          P_value = test$p.value
        )
      }, error = function(e) {
        tibble(
          DRUG_NAME = drug,
          PATHWAY_NAME = moa,
          refined_tumor_type = cell,
          Cluster = cluster,
          Correlation = NA_real_,
          P_value = NA_real_
        )
      })
    })
  }, .keep = TRUE) %>%
  bind_rows()


tmp <- filter(results, P_value <= 0.05, Correlation <0, DRUG_NAME %in% tc4_drugs, refined_tumor_type == "ESCA") ## Only one combination
plot_correlation(data, drug = "AZD2014", tumor_type = "ESCA", cluster = "3q25.2")

plot_correlation(data, drug = "AT13148", tumor_type = "ESCA", cluster = "3q26.33_markers")


tmp <- filter(results, P_value <= 0.05, Correlation <0, PATHWAY_NAME == "PI3K/MTOR signaling", refined_tumor_type == "ESCA") ## Two PI3Ki
plot_correlation(data, drug = "OSI-027", tumor_type = "ESCA", cluster = "3q27.2")

tmp <- filter(results, P_value <= 0.05, Correlation <0, refined_tumor_type == "ESCA")
table(tmp$PATHWAY_NAME, tmp$Cluster)


mat <- table(tmp$PATHWAY_NAME, tmp$Cluster)

# Dibujar el heatmap
library(pheatmap)

formatted_numbers <- matrix(as.character(round(mat, 0)), nrow = nrow(mat))

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = formatted_numbers,
         fontsize = 10,
         color = colorRampPalette(c("white", "#ae0000"))(100),
         breaks = seq(0, 4, length.out = 101),
         legend_breaks = c(0, 1, 2, 3, 4),
         legend_labels = c("0", "1", "2", "3", "4")
)

grid::grid.text("N drugs\nanticorrelated\nwith TC4 CNVs", x = 0.92, y = 0.94, gp = grid::gpar(fontsize = 9))


plot_correlation(data, drug = "OSI-027", tumor_type = "ESCA", cluster = "17q21.2_markers")


## Dotplot with 142 combinations of drug-primary tumor that show significant 
## anticorrelation of AUC with TC10 signature.
selected_drugs <- unique(tmp$DRUG_NAME)

d <- results %>%
  filter(DRUG_NAME %in% selected_drugs,
         refined_tumor_type == "ESCA") %>% 
  mutate(
    Significance = ifelse(P_value <= 0.05, "Significant", "Not Significant"),
    Corr_sign = ifelse(Correlation >= 0, "Positive", "Negative"),
    Corr_color = ifelse(Corr_sign == "Positive", "#db4646", "#0090ab")
  )


ggplot(d, aes(x = DRUG_NAME, y = Cluster)) +
  geom_point(aes(
    size = abs(Correlation),
    fill = Corr_color,
    color = Corr_color,
    shape = Significance
  ),
  stroke = 1.2  # para ver mejor el borde
  ) +
  scale_shape_manual(values = c("Significant" = 21, "Not Significant" = 1)) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_size(range = c(1, 7)) +
  facet_grid(. ~ PATHWAY_NAME, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 15),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "|Correlation|"
  )


ggsave("gdsc/gdsc_correlations_tc4.png", width = 20, height = 15, dpi = 300)

## Plot for figure
sigs_to_keep <- c("3q25.1", "3q25.2", "3q26.32", "3q26.33", 
                  "3q27.1", "3q27.2", "3q27.3", "3q28", "3q29")

tmp <- filter(results,
              Cluster %in% sigs_to_keep,
              P_value <= 0.05, Correlation <0, 
              refined_tumor_type == "ESCA")

selected_drugs <- unique(tmp$DRUG_NAME)



cluster_labels <- c(
  "3q25.1" = "3q25.1 (TM4SF1, PFN2)", 
  "3q25.2" = "3q25.2 (RAP2B)",
  "3q26.32" = "3q26.32", 
  "3q26.33" = "3q26.33 (ACTL6A)",
  "3q27.1" = "3q27.1 (PSMD2)", 
  "3q27.2" = "3q27.2 (MAP3K13)", 
  "3q27.3" = "3q27.3",
  "3q28" = "3q28 (CLDN1)", 
  "3q29" = "3q29 (HES1)"
)

d <- results %>%
  filter(Cluster %in% sigs_to_keep,
         DRUG_NAME %in% selected_drugs,
         refined_tumor_type == "ESCA") %>% 
  mutate(
    Significance = ifelse(P_value <= 0.05, "Significant", "Not Significant"),
    Corr_sign = ifelse(Correlation >= 0, "Positive", "Negative"),
    Corr_color = ifelse(Corr_sign == "Positive", "#db4646", "#0090ab")
  )


ggplot(d, aes(x = DRUG_NAME, y = Cluster)) +
  geom_point(aes(
    size = abs(Correlation),
    fill = Corr_color,
    color = Corr_color,
    shape = Significance
  ),
  stroke = 1.2  # para ver mejor el borde
  ) +
  scale_shape_manual(values = c("Significant" = 21, "Not Significant" = 1)) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_size(range = c(1, 7)) +
  facet_grid(. ~ PATHWAY_NAME, scales = "free_x", space = "free_x") +
  scale_y_discrete(labels = cluster_labels) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 15),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "|Correlation|"
  )


ggsave("gdsc/gdsc_correlations_tc4_figure.png", width = 20, height = 8, dpi = 300)

svg("gdsc/gdsc_correlations_tc4_figure.svg", width = 20, height = 8)
print(
  ggplot(d, aes(x = DRUG_NAME, y = Cluster)) +
    geom_point(aes(
      size = abs(Correlation),
      fill = Corr_color,
      color = Corr_color,
      shape = Significance
    ),
    stroke = 1.2  # para ver mejor el borde
    ) +
    scale_shape_manual(values = c("Significant" = 21, "Not Significant" = 1)) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_size(range = c(1, 7)) +
    facet_grid(. ~ PATHWAY_NAME, scales = "free_x", space = "free_x") +
    scale_y_discrete(labels = cluster_labels) + 
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 15),
      panel.grid = element_blank()
    ) +
    labs(
      x = NULL,
      y = NULL,
      size = "|Correlation|"
    )
)

dev.off() 

plot_correlation(data, drug = "AZD2014", tumor_type = "ESCA", cluster = "3q25.2")


## Plot for figure 2 (supplementary)
sigs_to_keep <- c("3q25.1", "3q25.2", "3q26.32", "3q26.33", 
                  "3q27.1", "3q27.2", "3q27.3", "3q28", "3q29")

tmp <- filter(results,
              Cluster %in% sigs_to_keep,
              P_value <= 0.05, Correlation >0, 
              refined_tumor_type == "ESCA")

selected_drugs <- unique(tmp$DRUG_NAME)



cluster_labels <- c(
  "3q25.1" = "3q25.1 (TM4SF1, PFN2)", 
  "3q25.2" = "3q25.2 (RAP2B)",
  "3q26.32" = "3q26.32", 
  "3q26.33" = "3q26.33 (ACTL6A)",
  "3q27.1" = "3q27.1 (PSMD2)", 
  "3q27.2" = "3q27.2 (MAP3K13)", 
  "3q27.3" = "3q27.3",
  "3q28" = "3q28 (CLDN1)", 
  "3q29" = "3q29 (HES1)"
)

d <- results %>%
  filter(Cluster %in% sigs_to_keep,
         DRUG_NAME %in% selected_drugs,
         refined_tumor_type == "ESCA") %>% 
  mutate(
    Significance = ifelse(P_value <= 0.05, "Significant", "Not Significant"),
    Corr_sign = ifelse(Correlation >= 0, "Positive", "Negative"),
    Corr_color = ifelse(Corr_sign == "Positive", "#db4646", "#0090ab")
  )


ggplot(d, aes(x = DRUG_NAME, y = Cluster)) +
  geom_point(aes(
    size = abs(Correlation),
    fill = Corr_color,
    color = Corr_color,
    shape = Significance
  ),
  stroke = 1.2  # para ver mejor el borde
  ) +
  scale_shape_manual(values = c("Significant" = 21, "Not Significant" = 1)) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_size(range = c(1, 7)) +
  facet_grid(. ~ PATHWAY_NAME, scales = "free_x", space = "free_x") +
  scale_y_discrete(labels = cluster_labels) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 15),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "|Correlation|"
  )


ggsave("gdsc/gdsc_correlations_tc4_figure2.png", width = 20, height = 8, dpi = 300)

## Plot for figure 1+2 (combined)
sigs_to_keep <- c("3q25.1", "3q25.2", "3q26.32", "3q26.33", 
                  "3q27.1", "3q27.2", "3q27.3", "3q28", "3q29")

tmp <- filter(results,
              Cluster %in% sigs_to_keep,
              P_value <= 0.05, 
              refined_tumor_type == "ESCA")

selected_drugs <- unique(tmp$DRUG_NAME)



cluster_labels <- c(
  "3q25.1" = "3q25.1 (TM4SF1, PFN2)", 
  "3q25.2" = "3q25.2 (RAP2B)",
  "3q26.32" = "3q26.32", 
  "3q26.33" = "3q26.33 (ACTL6A)",
  "3q27.1" = "3q27.1 (PSMD2)", 
  "3q27.2" = "3q27.2 (MAP3K13)", 
  "3q27.3" = "3q27.3",
  "3q28" = "3q28 (CLDN1)", 
  "3q29" = "3q29 (HES1)"
)

d <- results %>%
  filter(Cluster %in% sigs_to_keep,
         DRUG_NAME %in% selected_drugs,
         refined_tumor_type == "ESCA") %>% 
  mutate(
    Significance = ifelse(P_value <= 0.05, "Significant", "Not Significant"),
    Corr_sign = ifelse(Correlation >= 0, "Positive", "Negative"),
    Corr_color = ifelse(Corr_sign == "Positive", "#db4646", "#0090ab")
  )


ggplot(d, aes(x = DRUG_NAME, y = Cluster)) +
  geom_point(aes(
    size = abs(Correlation),
    fill = Corr_color,
    color = Corr_color,
    shape = Significance
  ),
  stroke = 1.2  # para ver mejor el borde
  ) +
  scale_shape_manual(values = c("Significant" = 21, "Not Significant" = 1)) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_size(range = c(1, 7)) +
  facet_grid(. ~ PATHWAY_NAME, scales = "free_x", space = "free_x") +
  scale_y_discrete(labels = cluster_labels) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 15),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "|Correlation|"
  )


ggsave("gdsc/gdsc_correlations_tc4_figure_combined.png", width = 26, height = 8, dpi = 300)

svg("gdsc/gdsc_correlations_tc4_figure_combined.svg", width = 26, height = 8)
  print(
    ggplot(d, aes(x = DRUG_NAME, y = Cluster)) +
      geom_point(aes(
        size = abs(Correlation),
        fill = Corr_color,
        color = Corr_color,
        shape = Significance
      ),
      stroke = 1.2  # para ver mejor el borde
      ) +
      scale_shape_manual(values = c("Significant" = 21, "Not Significant" = 1)) +
      scale_fill_identity() +
      scale_color_identity() +
      scale_size(range = c(1, 7)) +
      facet_grid(. ~ PATHWAY_NAME, scales = "free_x", space = "free_x") +
      scale_y_discrete(labels = cluster_labels) + 
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 15),
        panel.grid = element_blank()
      ) +
      labs(
        x = NULL,
        y = NULL,
        size = "|Correlation|"
      )
  )

dev.off() 
  
## Plot for everything, unfiltered
sigs_to_keep <- c("3q25.1", "3q25.2", "3q26.32", "3q26.33", 
                  "3q27.1", "3q27.2", "3q27.3", "3q28", "3q29")

tmp <- filter(results,
              P_value <= 0.05, 
              refined_tumor_type == "ESCA")

selected_drugs <- unique(tmp$DRUG_NAME)



cluster_labels <- c(
  "3q25.1" = "3q25.1 (TM4SF1, PFN2)", 
  "3q25.2" = "3q25.2 (RAP2B)",
  "3q26.32" = "3q26.32", 
  "3q26.33" = "3q26.33 (ACTL6A)",
  "3q27.1" = "3q27.1 (PSMD2)", 
  "3q27.2" = "3q27.2 (MAP3K13)", 
  "3q27.3" = "3q27.3",
  "3q28" = "3q28 (CLDN1)", 
  "3q29" = "3q29 (HES1)"
)

d <- results %>%
  filter(DRUG_NAME %in% selected_drugs,
         refined_tumor_type == "ESCA") %>% 
  mutate(
    Significance = ifelse(P_value <= 0.05, "Significant", "Not Significant"),
    Corr_sign = ifelse(Correlation >= 0, "Positive", "Negative"),
    Corr_color = ifelse(Corr_sign == "Positive", "#db4646", "#0090ab")
  )


ggplot(d, aes(x = DRUG_NAME, y = Cluster)) +
  geom_point(aes(
    size = abs(Correlation),
    fill = Corr_color,
    color = Corr_color,
    shape = Significance
  ),
  stroke = 1.2  # para ver mejor el borde
  ) +
  scale_shape_manual(values = c("Significant" = 21, "Not Significant" = 1)) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_size(range = c(1, 7)) +
  facet_grid(. ~ PATHWAY_NAME, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 15),
    panel.grid = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "|Correlation|"
  )


ggsave("gdsc/gdsc_correlations_tc4_figure_combined_unfiltered.png", width = 34, height = 16, dpi = 300)
ggsave("gdsc/gdsc_correlations_tc4_figure_combined_unfiltered.svg", width = 34, height = 16)
