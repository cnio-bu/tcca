library(stringr)
library(tidyverse)
library(edgeR)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(dplyr)


set.seed(1)
out.dir <- "/home/lmgonzalezb/Documents/bc-meta/SCellBow/TCGA/"
setwd(out.dir)


# --- Data ---
#TCGA expression data
expr.data <- data.frame(readRDS("pancancer_htseq_counts.rds"))
rownames(expr.data) <- expr.data$Ensembl_ID
expr.data$Ensembl_ID <- NULL

#TCGA metadata from Liu et al 2018
metadata <- read.table("Liu_2018.csv",
                       sep = ",", header = TRUE)
#Gene annotation
genes <- read.table("gencode.v22.annotation.gene.probeMap",
                    sep = '\t', header = TRUE)

# --- Code ---
# Select only primary tumor samples (code 01)
cols <- grep(pattern = "TCGA", colnames(expr.data), value = TRUE)
cols <- grep(pattern = "\\.01[A-Z]$", cols, value = TRUE)
expr.data <- expr.data[, cols]

# Edit column names so they match survival data
colnames(expr.data) <- str_replace_all(str_remove(colnames(expr.data),
                                                  pattern = "\\.01[A-Z]$"),
                                       pattern = "\\.", replacement = "-")

# Select common samples to both matrices
cols.to.keep <- intersect(colnames(expr.data), metadata$bcr_patient_barcode)
expr.data <- expr.data[, cols.to.keep]

# Transform log2(count + 1) to raw counts + round to the closest integer
expr.data <- round(2 ^ (expr.data) - 1, digits = 0)
expr.data[expr.data < 0] <- 0

#Generate a DGEList object
dge <- DGEList(counts = expr.data)

# Remove genes with counts consistently equal to zero or very low
keep <- filterByExpr(dge, min.prop = 0.6)
expr.data <- expr.data[keep, ]

# Scale normalization and voom
# dge <- calcNormFactors(dge, method = "TMM")
# dge.voom <- voom(dge, plot = TRUE)
# dge.logcpm <- dge.voom$E

# Annotate the genes
expr.data <- expr.data %>%
  as.data.frame() %>%
  rownames_to_column("id") %>% # Add a new column named "id" with the rownames
  left_join(genes[, c("id", "gene")], by = "id") %>% # Add gene name mapping column
  filter(!is.na(gene)) %>%
  dplyr::rename(Hugo_Symbol = gene, Entrez_Gene_Id = id) # Rename columns

# If there are duplicated genes, keep the max
expr.data <- aggregate(.~Hugo_Symbol, data = expr.data, max)

# Set gene symbol as rownames
rownames(expr.data) <- expr.data$Hugo_Symbol

# Remove id and symbol columns
gene_names <- rownames(expr.data)
expr.data <- expr.data[, -which(colnames(expr.data) %in%
                                            c("Hugo_Symbol", "Entrez_Gene_Id"))]

# Keep only GBM samples 
expr.data <- apply(expr.data, 2, as.numeric)

patients_gbm <- metadata %>%
  filter(type == "GBM" & bcr_patient_barcode %in% colnames(expr.data)) %>%
  pull(bcr_patient_barcode)

expr.data.gbm <- expr.data[, patients_gbm]
rownames(expr.data.gbm) <- gene_names

#Save expression data
write.table(expr.data.gbm, file = "survival_data_gbm.tsv", sep = "\t")


# Subset metadata from GBM samples
metadata.gbm <- metadata %>%
  filter(bcr_patient_barcode %in% patients_gbm) %>%
  mutate(subtype = type, 
         time = OS.time, 
         status = OS) %>%
  column_to_rownames(var = "bcr_patient_barcode") %>%
  dplyr::select(subtype, time, status) 

# Save metadata
write.table(metadata.gbm, file = "survival_metadata_gbm.tsv", sep = "\t")




