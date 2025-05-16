rm(list = ls())
library(stringr)
library(tidyverse)
library(edgeR)
library(GSVA)
library(GSEABase)
library(survival)
library(survminer)
library(survMisc)
library(glmnet)
library(ggpubr)

set.seed(1)
out.dir <-
  "/home/lmgonzalezb/Documents/bc-meta/sctherapy/marker_genes/"
setwd(out.dir)

# --- Functions ---
# Assign a tier to each sample according to its ssGSEA score
getHighLow <- function(x, top = 0.85, bottom = 0.15) {
  x_sorted <- sort(x, decreasing = TRUE)
  x_quantile <- quantile(x_sorted, c(top, bottom))
  high <- names(which(x_sorted > x_quantile[1]))
  low <- names(which(x_sorted < x_quantile[2]))
  return(list(high = high, low = low))
}

getResCox_output <- function(res.cox, variable) {
  res.cox <- summary(res.cox)
  p.value <- signif(res.cox$coef[variable, "Pr(>|z|)"], digits = 2)
  wald.statistic <- signif(res.cox$coef[variable, "z"], digits = 2)
  wald.test <- signif(as.numeric(res.cox$wald["test"]), digits = 2)
  p.value.test <-
    format(signif(as.numeric(res.cox$wald["pvalue"]), digits = 2), scientific = TRUE)
  beta <-
    signif(res.cox$coef[variable, 1], digits = 2)
  #coefficient beta
  HR <- signif(res.cox$coef[variable, 2], digits = 2)
  #exp(beta)
  HR.confint.lower <-
    signif(res.cox$conf.int[variable, "lower .95"], 2)
  HR.confint.upper <-
    signif(res.cox$conf.int[variable, "upper .95"], 2)
  HR <- paste0(HR, " (",
               HR.confint.lower, "-", HR.confint.upper, ")")
  
  res <- c(
    beta = beta,
    `HR (95% CI for HR)` = HR,
    wald.statistic = wald.statistic,
    p.value = p.value,
    wald.test = wald.test,
    p.value.test = p.value.test
  )
  
  return(res)
}

# --- Data ---
#TCGA expression data
expr.data <- data.frame(readRDS("./input_data_survival/pancancer_htseq_counts.rds"))
rownames(expr.data) <- expr.data$Ensembl_ID
expr.data$Ensembl_ID <- NULL

#TCGA metadata from Liu et al 2018
metadata <- read.table("./input_data_survival/Liu_2018.csv",
                       sep = ",",
                       header = TRUE)
#Gene annotation
genes <-
  read.table(
    "./input_data_survival/gencode.v22.annotation.gene.probeMap",
    sep = '\t',
    header = TRUE
  )

# --- Code ---
# Select only primary tumor samples (code 01)
cols <- grep(pattern = "TCGA", colnames(expr.data), value = TRUE)
cols <- grep(pattern = "\\.01[A-Z]$", cols, value = TRUE)
expr.data <- expr.data[, cols]

# Edit column names so they match survival data
colnames(expr.data) <-
  str_replace_all(
    str_remove(colnames(expr.data),
               pattern = "\\.01[A-Z]$"),
    pattern = "\\.",
    replacement = "-"
  )

# Select common samples to both matrices
cols.to.keep <-
  intersect(colnames(expr.data), metadata$bcr_patient_barcode)
expr.data <- expr.data[, cols.to.keep]

# Transform log2(count + 1) to raw counts + round to the closest integer
expr.data <- round(2 ^ (expr.data) - 1, digits = 0)
expr.data[expr.data < 0] <- 0

#Generate a DGEList object
dge <- DGEList(counts = expr.data)

# Remove genes with counts consistently equal to zero or very low
keep <- filterByExpr(dge, min.prop = 0.6)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Scale normalization and voom
dge <- calcNormFactors(dge, method = "TMM")
dge.voom <- voom(dge, plot = TRUE)
dge.logcpm <- dge.voom$E

# Annotate the genes
expr.data <- dge.logcpm %>%
  as.data.frame() %>%
  rownames_to_column("id") %>% # Add a new column named "id" with the rownames
  left_join(genes[, c("id", "gene")], by = "id") %>% # Add gene name mapping column
  filter(!is.na(gene)) %>%
  dplyr::rename(Hugo_Symbol = gene, Entrez_Gene_Id = id) # Rename columns

# If there are duplicated genes, keep the max
expr.data <- aggregate(. ~ Hugo_Symbol, data = expr.data, max)

# Set gene symbol as rownames
rownames(expr.data) <- expr.data$Hugo_Symbol

# Remove id and symbol columns
expr.data <- as.matrix(expr.data[,-which(colnames(expr.data) %in%
                                           c("Hugo_Symbol", "Entrez_Gene_Id"))])

# Set to numeric
expr.data <-
  matrix(
    as.numeric(expr.data),
    ncol = ncol(expr.data),
    byrow = FALSE,
    dimnames = list(rownames(expr.data), colnames(expr.data))
  )

#Save
saveRDS(expr.data, file = "expr.data_tcga.rds")

# Create list with the gene set
collection <- getGmt("marker_sigs_clusters.gmt",
                     geneIdType = GenenameIdentifier(),
                     collectionType = ComputedCollection())

# Compute GSVA for the genesets (this step required high RAM resources so it was
# run in the CNIO HPC cluster)
gsva <- gsva(
  expr = expr.data,
  gset.idx.list = collection,
  method = "ssgsea",
  ssgsea.norm = TRUE,
  min.sz = 10,
  max.sz = 1800,
  parallel.sz = 60
)

saveRDS(gsva, "gsva.rds")



# Create a prognostic model of survival based on the functional gene expression
# signatures:
gsva <- readRDS("gsva.rds")
# We join the table of unscaled gsva scores to the metadata information.
gsva <- as.data.frame(t(gsva))
signatures <- colnames(gsva)
gsva["bcr_patient_barcode"] <- rownames(gsva)
gsva_metadata <-
  left_join(gsva, metadata, by = "bcr_patient_barcode")


# Compute Cox Proportional-Hazards Model per gene set (adjusting by cancer type,
#gender, age and tumor stage).
fix_cov <-
  c("type",
    "gender",
    "age_at_initial_pathologic_diagnosis",
    "clinical_stage")

# Cox regression analysis with continuous variable (GSVA score per signature)
# using PFI (Progression Free Survival) as response variable.
formulas_PFI <- sapply(signatures,
                       function(x)
                         as.formula(paste(
                           'Surv(PFI.time, PFI)~',
                           paste(c(fix_cov, x), collapse = " + ")
                         )))

models_PFI <-
  lapply(formulas_PFI, function(x) {
    coxph(x, data = gsva_metadata)
  })

results_PFI <-
  as.data.frame(t(mapply(getResCox_output, models_PFI, signatures)))
saveRDS(models_PFI, "cox_models_pfi.rds")
write.table(
  results_PFI,
  file = "cox_pfi.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)


# Cox regression analysis with continuous variable (GSVA score per signature)
# using OS (Overall Survival) as response variable.
formulas_OS <- sapply(signatures,
                      function(x)
                        as.formula(paste(
                          'Surv(OS.time, OS)~',
                          paste(c(fix_cov, x), collapse = " + ")
                        )))

models_OS <-
  lapply(formulas_OS, function(x) {
    coxph(x, data = gsva_metadata)
  })

results_OS <-
  as.data.frame(t(mapply(getResCox_output, models_OS, signatures)))
saveRDS(models_OS, "cox_models_os.rds")
write.table(
  results_OS,
  file = "cox_os.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)




### Perform the analysis only with cluster 4 and ESCA cancer
signatures <-
  grep("Cluster04",
       colnames(gsva),
       ignore.case = TRUE,
       value = TRUE)
gsva_metadata_clust4 <- gsva_metadata %>%
  filter(type == "ESCA")


# Compute Cox Proportional-Hazards Model per gene set (adjusting by cancer type,
#gender, age and tumor stage).
fix_cov <- c("gender",
             "age_at_initial_pathologic_diagnosis",
             "clinical_stage")

# Cox regression analysis with continuous variable (GSVA score per signature)
# using PFI (Progression Free Survival) as response variable.
formulas_PFI <- sapply(signatures,
                       function(x)
                         as.formula(paste(
                           'Surv(PFI.time, PFI)~',
                           paste(c(fix_cov, x), collapse = " + ")
                         )))

models_PFI <-
  lapply(formulas_PFI, function(x) {
    coxph(x, data = gsva_metadata_clust4)
  })

results_PFI <-
  as.data.frame(t(mapply(getResCox_output, models_PFI, signatures)))
saveRDS(models_PFI, "cox_models_pfi_clust4.rds")
write.table(
  results_PFI,
  file = "cox_pfi_clust4.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)


# Cox regression analysis with continuous variable (GSVA score per signature)
# using OS (Overall Survival) as response variable.
formulas_OS <- sapply(signatures,
                      function(x)
                        as.formula(paste(
                          'Surv(OS.time, OS)~',
                          paste(c(fix_cov, x), collapse = " + ")
                        )))

models_OS <-
  lapply(formulas_OS, function(x) {
    coxph(x, data = gsva_metadata_clust4)
  })

results_OS <-
  as.data.frame(t(mapply(getResCox_output, models_OS, signatures)))
saveRDS(models_OS, "cox_models_os_clust4.rds")
write.table(
  results_OS,
  file = "cox_os_clust4.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)




### Perform the analysis only with cluster 5 and GBM cancer
signatures <-
  grep("Cluster05",
       colnames(gsva),
       ignore.case = TRUE,
       value = TRUE)
gsva_metadata_clust5 <- gsva_metadata %>%
  filter(type == "GBM")


# Compute Cox Proportional-Hazards Model per gene set (adjusting by cancer type,
#gender, age and tumor stage).
fix_cov <- c("gender", "age_at_initial_pathologic_diagnosis")

# Cox regression analysis with continuous variable (GSVA score per signature)
# using PFI (Progression Free Survival) as response variable.
formulas_PFI <- sapply(signatures,
                       function(x)
                         as.formula(paste(
                           'Surv(PFI.time, PFI)~',
                           paste(c(fix_cov, x), collapse = " + ")
                         )))

models_PFI <-
  lapply(formulas_PFI, function(x) {
    coxph(x, data = gsva_metadata_clust5)
  })

results_PFI <-
  as.data.frame(t(mapply(getResCox_output, models_PFI, signatures)))
saveRDS(models_PFI, "cox_models_pfi_clust5.rds")
write.table(
  results_PFI,
  file = "cox_pfi_clust5.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)


# Cox regression analysis with continuous variable (GSVA score per signature)
# using OS (Overall Survival) as response variable.
formulas_OS <- sapply(signatures,
                      function(x)
                        as.formula(paste(
                          'Surv(OS.time, OS)~',
                          paste(c(fix_cov, x), collapse = " + ")
                        )))

models_OS <-
  lapply(formulas_OS, function(x) {
    coxph(x, data = gsva_metadata_clust5)
  })

results_OS <-
  as.data.frame(t(mapply(getResCox_output, models_OS, signatures)))
saveRDS(models_OS, "cox_models_os_clust5.rds")
write.table(
  results_OS,
  file = "cox_os_clust5.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)


# Plot KM curve for cluster 10
res.cut <- surv_cutpoint(gsva_metadata,
                         time = "PFI.time", event = "PFI",
                         variables = "Cluster10_top100_UP")
gsva_metadata$clust10_pfi <- ifelse(gsva_metadata$Cluster10_top100_UP < res.cut$cutpoint$cutpoint, "Low", "High")
gsva_metadata$clust10_pfi <- factor(gsva_metadata$clust10_pfi)

# Set the reference levels for each categorical variable. We try to use as reference 
# levels the ones with lower malignancy.
gsva_metadata$clust10_pfi<- relevel(gsva_metadata$clust10_pfi, ref = "Low")

ggsurv1 <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ clust10_pfi, data = gsva_metadata),
                      conf.int = T, risk.table = TRUE,
                      palette = c("#28cfb8", "#FF6663"),
                      legend.labs = c("Low score", "High score"))

png("clust10_pfi_curve.png", width = 8, height = 7, units = "in", res = 500)
ggsurv1
dev.off()

# The same for OS
res.cut <- surv_cutpoint(gsva_metadata,
                         time = "OS.time", event = "OS",
                         variables = "Cluster10_top100_UP")
gsva_metadata$clust10_os <- ifelse(gsva_metadata$Cluster10_top100_UP < res.cut$cutpoint$cutpoint, "Low", "High")
gsva_metadata$clust10_os <- factor(gsva_metadata$clust10_os)

gsva_metadata$clust10_os <- relevel(gsva_metadata$clust10_os, ref = "Low")

ggsurv1 <- ggsurvplot(survfit(Surv(OS.time, OS) ~ clust10_os, data = gsva_metadata),
                      conf.int = T, risk.table = TRUE,
                      palette = c("#28cfb8", "#FF6663"),
                      legend.labs = c("Low score", "High score"))

png("clust10_os_curve.png", width = 8, height = 7, units = "in", res = 500)
ggsurv1
dev.off()




## Check survival for amplified marker genes in cluster 4 (ESCA)
genes_band1 <- c("S100A14", "S100A16", "S100A2", "S100A7", "S100A8", "SPRR1B")
genes_band2 <- c("JUP", "KRT14", "KRT15", "KRT16", "KRT17", "TNS4")

gsva <- gsva(
  expr = expr.data,
  gset.idx.list = list(genes_band1_c4 = genes_band1, genes_band2_c4 = genes_band2),
  method = "ssgsea",
  ssgsea.norm = TRUE,
  min.sz = 3,
  max.sz = 1800,
  parallel.sz = 2
)

gsva <- as.data.frame(t(gsva))
signatures <- colnames(gsva)
gsva["bcr_patient_barcode"] <- rownames(gsva)
gsva_metadata <-
  left_join(gsva, metadata, by = "bcr_patient_barcode")
gsva_metadata_clust4 <- gsva_metadata %>%
  filter(type == "ESCA")

# Compute survival
fix_cov <- c("gender","age_at_initial_pathologic_diagnosis", "clinical_stage")

# Cox regression analysis with continuous variable (GSVA score per signature)
# using PFI (Progression Free Survival) as response variable.
formulas_PFI <- sapply(signatures,
                       function(x)
                         as.formula(paste(
                           'Surv(PFI.time, PFI)~',
                           paste(c(fix_cov, x), collapse = " + ")
                         )))

models_PFI <-
  lapply(formulas_PFI, function(x) {
    coxph(x, data = gsva_metadata_clust4)
  })


# The same using OS
formulas_OS <- sapply(signatures,
                      function(x)
                        as.formula(paste(
                          'Surv(OS.time, OS)~',
                          paste(c(fix_cov, x), collapse = " + ")
                        )))

models_OS <-
  lapply(formulas_OS, function(x) {
    coxph(x, data = gsva_metadata_clust4)
  })






## Find the optimal gene signature for cluster 10 starting from 25 genes
gene_sets <- lapply(original_signature, function(gene) setdiff(original_signature, gene))
names(gene_sets) <- paste0("no_", original_signature)

# Also include the full signature if you want to compare
gene_sets$full_signature <- original_signature

# Run GSVA (or ssGSEA) on all signatures in one shot
gsva_scores <- gsva(expression_matrix, gene_sets, method = "ssgsea", verbose = FALSE)

library(survival)

results <- data.frame(
  signature = character(),
  HR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (sig_name in rownames(gsva_scores)) {
  gsva_metadata$gsva_score <- as.vector(gsva_scores[sig_name, ])
  
  model <- coxph(Surv(PFI.time, PFI) ~ gsva_score + gender + age + clinical_stage, data = gsva_metadata)
  s <- summary(model)
  
  results <- rbind(results, data.frame(
    signature = sig_name,
    HR = s$coefficients["gsva_score", "exp(coef)"],
    CI_lower = s$conf.int["gsva_score", "lower .95"],
    CI_upper = s$conf.int["gsva_score", "upper .95"],
    pvalue = s$coefficients["gsva_score", "Pr(>|z|)"]
  ))
}
