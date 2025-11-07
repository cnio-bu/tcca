library(Seurat)
library(BPCells)
library(UCell)
library(dplyr)
library(tidyverse)
library(ggplot2)

setwd("/storage/scratch01/shared/projects/bc-meta/single_cell/")

##---------------------------------Functions ---------------------------------##
read.gmt <- function(gmt_file) {
  sigs_list <- list()
  sigs <- scan(gmt_file, what = character(), sep = "\n")
  for (sig in sigs) {
    sig <- unlist(strsplit(sig, "\t"))
    sig <- unique(sig[nzchar(sig)])
    sigs_list[[sig[1]]] <- sig[3:length(sig)]
  }
  return(sigs_list)
}

# Read the gmt file with the VAR score gene sets.
gmt_var <- read.gmt(
  "/home/mgonzalezb/bc-meta/var_score/gene_sets/selected_functional_signatures.gmt"
)

# Join bidirectional signatures in the same vector adding '+' to the genes in the 
# upregulared gene set and '-' to the genes in the downregulated gene set.
sig_names <- sub("_(UP|DOWN|DN)$", "", names(gmt_var))
bisigs <- sig_names[duplicated(sig_names)]

for (sig in unique(sig_names)) {
  up <- paste0(sig, "_UP")
  down <-  if (paste0(sig, "_DOWN") %in% names(gmt_var))
    paste0(sig, "_DOWN")
  else
    paste0(sig, "_DN")
  
  if (sig %in% bisigs) {
    gmt_var[[up]] <- paste0(gmt_var[[up]], "+")
    gmt_var[[down]] <- paste0(gmt_var[[down]], "-")
    gmt_var[[sig]] <- c(gmt_var[[up]], gmt_var[[down]])
    gmt_var[[up]] <- NULL
    gmt_var[[down]] <- NULL
    
  } else{
    if (up %in% names(gmt_var)) {
      gmt_var[[up]] <- paste0(gmt_var[[up]], "+")
    }
    
    if (down %in% names(gmt_var)) {
      gmt_var[[down]] <- paste0(gmt_var[[down]], "-")
    }
  }
}


# Load level 2 Seurat object
seu <- readRDS(
  "seurat/v5/lvl2/seu_lvl2.rds"
)
malignant <- subset(seu, subset = malignancy == TRUE)
colnames(malignant) <- paste0("c", c(1:ncol(malignant)))

# Compute UCell score for the VAR signatures
seu_ucell <- AddModuleScore_UCell(malignant, features = gmt_var)

# Save UCell scores separately as a matrix
sig_names <- paste0(names(gmt_var), "_UCell")
scores_ucell <- seu_ucell@meta.data[, sig_names]
write.table(scores_ucell, "varscore/ucell_varsigs.tsv")

# Compute A and VR scores per cell
a.model <- readRDS("/home/mgonzalezb/bc-meta/var_score/multivariable_Res.cox_A_score.rds")
vr.model <- readRDS("/home/mgonzalezb/bc-meta/var_score/multivariable_Res.cox_VR_score.rds")
colnames(scores_ucell) <- gsub("_UCell", "", colnames(scores_ucell))
ascore_raw <- predict(a.model, newdata = scores_ucell, type = "lp")
vrscore_raw <- predict(vr.model, newdata = scores_ucell, type = "lp")
malignant$study_sample <- paste0(malignant$study, "_", malignant$sample)

var_scores <- data.frame(A_score_raw = ascore_raw, VR_score_raw = vrscore_raw)
identical(rownames(var_scores), rownames(malignant@meta.data))
var_scores <- cbind(malignant@meta.data[, "study_sample", drop = FALSE], var_scores)
head(var_scores)

# Sample-wise scale A and VR scores
var_scores <- var_scores %>%
  rownames_to_column("RowName") %>%
  group_by(study_sample) %>%
  mutate(A_zscore = (A_score_raw - mean(A_score_raw))/sd(A_score_raw),
         VR_zscore = (VR_score_raw - mean(VR_score_raw))/sd(VR_score_raw),
         A_scaled = A_zscore/max(abs(A_zscore)),
         VR_scaled = VR_zscore/max(abs(VR_zscore)),
         VARscore_raw_subs = A_zscore - VR_zscore,
         VARscore_scaled_subs = (VARscore_raw_subs - min(VARscore_raw_subs))/
         (max(VARscore_raw_subs) - min(VARscore_raw_subs)),
         VARscore_raw_sum = A_zscore + VR_zscore,
         VARscore_scaled_sum = (VARscore_raw_sum - min(VARscore_raw_sum))/
         (max(VARscore_raw_sum) - min(VARscore_raw_sum))) %>%
  ungroup() %>%
  column_to_rownames("RowName") %>%
  mutate(across(!study_sample, ~ signif(.x, 3))) %>%
  as.data.frame()

write.table(var_scores, "varscore/var_scores.tsv")

