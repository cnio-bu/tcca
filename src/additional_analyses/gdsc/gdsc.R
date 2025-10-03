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
  
  test <- cor.test(df_plot[[cluster]], df_plot$AUC, method = "pearson")
  corr_val <- round(test$estimate, 3)
  p_val <- test$p.value
  p_val_formatted <- ifelse(p_val < 0.001, "< 0.001", paste0("= ", round(p_val, 4)))
  
  if (unique(tumor_type == "all")) {
    ggplot(df_plot, aes_string(x = cluster, y = "AUC", color = "refined_tumor_type")) +
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
    
    ggplot(df_plot, aes_string(x = cluster, y = "AUC", color = "refined_tumor_type")) +
      geom_point() +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      labs(
        title = paste0("Pearson correlation: ", cluster, " and AUC"),
        subtitle = paste0("Drug: ", drug, " | ", tumor_text, "\n", "Pearson r = ", corr_val, ", p ", p_val_formatted),
        x = cluster,
        y = "AUC"
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

markers <- read_gmt_list("genes_by_cluster.gmt")

## Adapt signatures to keep only genes expressed in the dataset

expr_matrix <- GetAssayData(seu, slot = "counts")  # o "data" para valores normalizados
expr_genes <- rownames(expr_matrix)[Matrix::rowSums(expr_matrix > 0) > 0]

markers <- lapply(markers, function(genes_vector) {
  intersect(genes_vector, expr_genes)
})


seu <- AddModuleScore_UCell(seu, features=markers, name = NULL)

seu <- AddModuleScore_UCell(seu, features=markers, name = NULL)

## Check undetected genes
undetected_genes <- c(
  "CD24", "PLPP2", "PRELID3B", "RACK1", "NOP53", "ATP5F1E", "ATP5MG",
  "SELENOK", "SELENOS", "MESD", "REX1BD", "SEM1", "RFLNB", "VSIR",
  "FYB1", "GRK2", "ADGRE5", "RAB5IF", "STMP1", "ATP5MC2", "RTRAF"
)

## Undetected genes are mainly (13/21) from TC8
length(intersect(markers$Cluster10_UP, undetected_genes))

## Get dataframe of averages

tc_sigs <- seu@meta.data %>%
  dplyr::select(CELL_LINE_NAME, refined_tumor_type, starts_with("Cluster")) %>%
  group_by(CELL_LINE_NAME, refined_tumor_type) %>%
  summarise(across(starts_with("Cluster"), mean, na.rm = TRUE), .groups = "drop")

data <- left_join(tc_sigs, gdsc)


tc10_drugs <- c(
  "AZD8055",
  "Bortezomib",
  "Carfilzomib",
  "Danusertib",
  "Daunorubicin",
  "Delanzomib",
  "Elesclomol",
  "Epirubicin",
  "Ixazomib",
  "JNJ-26481585",
  "Mitoxantrone",
  "MLN0128",
  "NVP-bez235",
  "NVP-tae226",
  "Panobinostat",
  "Romidepsin",
  "SN-38"
)

tc10_drugs <- intersect(tc10_drugs, data$DRUG_NAME)

tc10_tumors <- c(
  "HNSC",  # head and neck
  "COAD", "READ",  # colon / colorectal
  "SKCM",  # skin
  "ESCA",  # esophageal
  "BRCA",  # breast
  "LUAD", "LUSC", "LCLC", "SCLC",  # lung
  "PRAD",  # prostate
  "OV",    # ovarian
  "PAAD",  # pancreas
  "GBM", "LGG"  # brain
)

## Calculate correlations with all tumor types
cluster_cols <- grep("Cluster\\d+_UP", colnames(data), value = TRUE)

results <- data %>%
  group_by(DRUG_NAME) %>%
  group_map(~ {
    drug <- unique(.x$DRUG_NAME)
    moa <- unique(.x$PATHWAY_NAME)
    map_dfr(cluster_cols, function(cluster) {
      test <- cor.test(.x[[cluster]], .x$AUC, method = "pearson")
      tibble(
        DRUG_NAME = drug,
        PATHWAY_NAME = moa,
        Cluster = cluster,
        Correlation = test$estimate,
        P_value = test$p.value
      )
    })
  }, .keep = TRUE) %>%
  bind_rows()

tmp <- filter(results, P_value < 0.05, Cluster == "Cluster10_UP") # + DRUG_NAME %in% tc10_drugs returns 0
plot_correlation(data, drug = "5-azacytidine", tumor_type = "all", cluster = "Cluster10_UP")


## Calculate correlations with main tumor types from cluster 10
cluster_cols <- grep("Cluster\\d+_UP", colnames(data), value = TRUE)


results2 <- data %>%
  filter(refined_tumor_type %in% tc10_tumors) %>%
  group_by(DRUG_NAME) %>%
  group_map(~ {
    drug <- unique(.x$DRUG_NAME)
    moa <- unique(.x$PATHWAY_NAME)
    map_dfr(cluster_cols, function(cluster) {
      test <- cor.test(.x[[cluster]], .x$AUC, method = "pearson")
      tibble(
        DRUG_NAME = drug,
        PATHWAY_NAME =moa,
        Cluster = cluster,
        Correlation = test$estimate,
        P_value = test$p.value
      )
    })
  }, .keep = TRUE) %>%
  bind_rows()

tmp <- filter(results2, P_value < 0.05, Cluster == "Cluster10_UP") # + DRUG_NAME %in% tc10_drugs returns 0
plot_correlation(data, drug = "5-azacytidine", tumor_type = tc10_tumors, cluster = "Cluster10_UP")



## Por tumor type

# Calcular correlación por combinación DRUG_NAME + CELL_LINE_NAME
results3 <- data %>%
  group_by(DRUG_NAME, refined_tumor_type) %>%
  group_map(~ {
    drug <- unique(.x$DRUG_NAME)
    cell <- unique(.x$refined_tumor_type)
    moa <- unique(.x$PATHWAY_NAME)
    
    # Verificamos que haya suficiente variabilidad y datos
    if (nrow(.x) < 3) return(NULL)
    
    map_dfr(cluster_cols, function(cluster) {
      # Control de errores por si hay NA o varianza 0
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


tmp <- filter(results3, Cluster == "Cluster10_UP", P_value <= 0.05, Correlation <0) #, DRUG_NAME %in% tc10_drugs) 
table(tmp$PATHWAY_NAME)

tmp <- filter(results3, Cluster == "Cluster04_UP", P_value <= 0.05, Correlation <0, refined_tumor_type == "ESCA") ## 0 OBSERVATIONS wtfff 
table(tmp$PATHWAY_NAME)

tmp <- filter(results3, Cluster == "Cluster05_UP", P_value <= 0.05, Correlation <0) 
table(tmp$PATHWAY_NAME)

tmp <- filter(results3, Cluster == "Cluster06_UP", P_value <= 0.05, Correlation <0) 
table(tmp$PATHWAY_NAME)

tmp <- filter(results3, P_value <= 0.05, Correlation <0, DRUG_NAME %in% all_drugs) 

tmp <- filter(results3, P_value <= 0.05, Correlation <0, DRUG_NAME %in% cluster_10, Cluster == "Cluster10_UP") 
plot_correlation(data, drug = "Bortezomib", tumor_type = "SARC", cluster = "Cluster02_UP")

tmp <- filter(results3, P_value <= 0.05, Correlation <0, Cluster == "Cluster10_UP") 


plot_correlation(data, drug = "Alpelisib", tumor_type = "COAD", cluster = "Cluster10_UP")
plot_correlation(data, drug = "AZD8186", tumor_type = "ESCA", cluster = "Cluster04_UP")
plot_correlation(data, drug = "BI-2536", tumor_type = "GBM", cluster = "Cluster05_UP")

tmp <- filter(results3, refined_tumor_type == "BRCA", DRUG_NAME == "Epirubicin") 
tmp <- filter(results3, refined_tumor_type == "GBM", DRUG_NAME == "AZD8055") 
tmp <- filter(results3, refined_tumor_type == "LUSC", DRUG_NAME == "Mitoxantrone") 

plot_correlation(data, drug = "Epirubicin", tumor_type = "BRCA", cluster = "Cluster10_UP")

svg("epirubicin.svg", width=8, height=6)

plot_correlation(data, drug = "Epirubicin", tumor_type = "BRCA", cluster = "Cluster10_UP")

dev.off()




plot_correlation(data, drug = "AZD8055", tumor_type = "GBM", cluster = "Cluster10_UP")
plot_correlation(data, drug = "Mitoxantrone", tumor_type = "LUSC", cluster = "Cluster10_UP")

plot_correlation(data, drug = "Veliparib", tumor_type = "COAD", cluster = "Cluster10_UP")



tc10_sig <- filter(results3, P_value <= 0.05, Correlation <0, Cluster == "Cluster10_UP") %>% mutate(combo = paste0(DRUG_NAME, "_", refined_tumor_type)) %>% dplyr::select(combo) %>% unique() %>% unlist()
other_ns <- filter(results3, P_value > 0.05, Cluster != "Cluster10_UP") %>% mutate(combo = paste0(DRUG_NAME, "_", refined_tumor_type)) %>% dplyr::select(combo) %>% unique() %>% unlist()

drugs_of_interest <- intersect(tc10_sig, other_ns)

doi <- results3 %>%
  mutate(combo = paste0(DRUG_NAME, "_", refined_tumor_type)) %>%
  filter(combo %in% drugs_of_interest)


tmp <- filter(doi, Cluster == "Cluster10_UP")

mat <- table(tmp$PATHWAY_NAME, tmp$refined_tumor_type)

# Dibujar el heatmap
library(pheatmap)

formatted_numbers <- matrix(as.character(round(mat, 0)), nrow = nrow(mat))

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = formatted_numbers,
         fontsize = 10,
         color = colorRampPalette(c("white", "#ae0000"))(100),
         breaks = seq(0, 5, length.out = 101),
         legend_breaks = c(0, 1, 2, 3, 4, 5),
         legend_labels = c("0", "1", "2", "3", "4", "5+")
         )

grid::grid.text("N drugs\nanticorrelated\nwith TC10", x = 0.94, y = 0.94, gp = grid::gpar(fontsize = 9))

## Dotplot with 142 combinations of drug-primary tumor that show significant 
## anticorrelation of AUC with TC10 signature.

selected_drugs <- unique(doi$DRUG_NAME)
selected_tt <- unique(doi$refined_tumor_type)

d <- results3 %>%
  filter(Cluster == "Cluster10_UP", 
         DRUG_NAME %in% selected_drugs) %>% 
         #refined_tumor_type %in% selected_tt) %>%
  mutate(
    Significance = ifelse(P_value <= 0.05, "Significant", "Not Significant"),
    Corr_sign = ifelse(Correlation >= 0, "Positive", "Negative"),
    Corr_color = ifelse(Corr_sign == "Positive", "#db4646", "#0090ab")
  )


ggplot(d, aes(x = DRUG_NAME, y = refined_tumor_type)) +
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


ggsave("gdsc/gdsc_correlations_142.png", width = 30, height = 11, dpi = 300)

svg("gdsc/gdsc_correlations_142.svg", width = 30, height = 11)
print(
  ggplot(d, aes(x = DRUG_NAME, y = refined_tumor_type)) +
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
)

dev.off() 

## Keep only high TC10 samples

hightc10 <- quantile(data$Cluster10_UP, 0.25)

## Calculate correlations with all tumor types
cluster_cols <- grep("Cluster\\d+_UP", colnames(data), value = TRUE)

results4 <- data %>%
  filter(Cluster10_UP <= hightc10) %>%
  group_by(DRUG_NAME) %>%
  group_map(~ {
    drug <- unique(.x$DRUG_NAME)
    moa <- unique(.x$PATHWAY_NAME)
    map_dfr(cluster_cols, function(cluster) {
      test <- cor.test(.x[[cluster]], .x$AUC, method = "pearson")
      tibble(
        DRUG_NAME = drug,
        PATHWAY_NAME = moa,
        Cluster = cluster,
        Correlation = test$estimate,
        P_value = test$p.value
      )
    })
  }, .keep = TRUE) %>%
  bind_rows()

tmp <- filter(results4, P_value < 0.05, Cluster == "Cluster10_UP", DRUG_NAME %in% tc10_drugs)
plot_correlation(data %>% filter(Cluster10_UP <= hightc10), drug = "Epirubicin", tumor_type = "all", cluster = "Cluster10_UP")
