library(clusterProfiler)
library(org.Hs.eg.db)
library(BPCells)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(AnnotationHub)
library(dplyr)
library(purrr)
library(broom)
library(stringr)
library(enrichplot)
library(grid)


setwd("/local/bc_meta/")

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

rm(seu, expr_matrix, gdsc, metadata)
gc()

## GO Enrichment

markers_go <- markers[1:12] ## Cross with TC4 biomarkers are too few genes per signature

markers_go <- lapply(markers_go, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
})

go_results <- lapply(markers_go, function(genes) {
  enrichGO(gene = genes,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable = TRUE)
})


for (i in 1:length(go_results)) {
  print(
    dotplot(go_results[[i]], showCategory = 10) + 
      ggtitle(names(go_results)[i])
  )
}

dotplot(go_results[["3q25.2"]], showCategory = 100) + 
  ggtitle("3q25.2")

go3q25 <- as.data.frame(go_results[["3q25.2"]])
go3q25 <- go3q25[grepl("morpho|epith|emt|EMT|transition|mesen|stress", go3q25$Description, ignore.case = TRUE), ]

dotplot(go_results[["3q25.2"]], showCategory = 100) + 
  ggtitle("3q25.2")

library(dplyr)

go_df <- bind_rows(lapply(names(go_results), function(name) {
  res <- go_results[[name]]
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  
  as.data.frame(res) %>%
    mutate(GeneSet = name)
}))

go_df <- go_df %>%
  mutate(
    log10padj = -log10(p.adjust),
    GeneRatioNum = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
  )

top_terms <- go_df %>%
  group_by(GeneSet) %>%
  arrange(p.adjust, desc(GeneRatioNum), .by_group = TRUE) %>%  # ordenar por p.adjust y luego por enriquecimiento
  mutate(rank = row_number()) %>%
  filter(rank <= 30) %>%
  ungroup()

top_terms <- filter(go_df, ID %in% top_terms$ID)

top_terms_ordered <- top_terms %>%
  group_by(GeneSet) %>%
  arrange(desc(log10padj), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(
    GeneSet = factor(GeneSet, levels = unique(GeneSet)),
    Description = factor(Description, levels = rev(unique(Description)))  # Para que vaya de izquierda a derecha
  )

ggplot(top_terms_ordered, aes(x = Description, y = GeneSet)) +
  geom_point(aes(size = log10padj, color = GeneRatioNum)) +
  scale_color_viridis_c(option = "plasma", name = "GeneRatio") +
  scale_size(name = "-log10(p.adjust)") +
  theme_bw() +
  labs(
    x = "GO Term",
    y = "Gene Set",
    title = "GO Enrichment Dotplot"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    plot.margin = unit(c(1, 1, 1, 15), "lines")
  )

ggsave("gdsc/gdsc_goenrichment_tc4.png", width = 30, height = 7, dpi = 200)