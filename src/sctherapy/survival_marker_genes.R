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
#library(survcomp)
library(forestmodel)
library(openxlsx)
library(car)

set.seed(1)
out.dir <-
  "/home/lmgonzalezb/Documents/bc-meta/sctherapy/marker_genes/"
setwd(out.dir)

# --- Functions ---
# Assign a tier to each sample according to its ssGSEA score
getHighLow <- function(x, top = 0.75, bottom = 0.25) {
  x_sorted <- sort(x, decreasing = TRUE)
  x_quantile <- quantile(x_sorted, c(top, bottom))
  high <- names(which(x_sorted > x_quantile[1]))
  low <- names(which(x_sorted < x_quantile[2]))
  return(list(high = high, low = low))
}

# Other way of selecting the optimal cutpoint according to its ssGSEA score
opt_cutpoint <- function(data, surv_param, continuous_var){
  surv_event <- surv_param
  surv_time <- paste0(surv_event, ".time")
  
  cut_values <- quantile(data[[continuous_var]], probs = seq(0.2, 0.8, 0.01))
  results <- data.frame(cutpoint = cut_values, pval = NA, cindex = NA)
  
  for (i in seq_along(cut_values)){
    cut <- cut_values[i]
    data$group <- ifelse(data[[continuous_var]] > cut, "High", "Low")
    data$group <- factor(data$group, levels = c("Low", "High"))
    
    surv_obj <- Surv(data[[surv_time]], data[[surv_event]])
    
    # Log-rank test for KM curves
    pval <- survdiff(surv_obj ~ group, data = data)
    pval <- pval$pvalue
    
    # C-index from Cox model
    cox_model <- coxph(surv_obj ~ group, data = data)
    c_index <- summary(cox_model)$concordance[1]
    
    results$pval[i] <- pval
    results$cindex[i] <- c_index
  }
  opt_cutpoint <- results[which.min(results$pval/results$cindex), ]
  return(list(optimal_cutpoint = opt_cutpoint, results = results))
}

# Plot C-index and p-value against possible cutpoints
opt_cutpoint_plot <- function(opt_cutpoint_output) {
  # Extract results and calculate -log10(p)
  cutpoint_values <- opt_cutpoint_output$results
  cutpoint_values$neg_log10_p <- -log10(cutpoint_values$pval)
  
  # Create the plot
  plot <- ggplot(cutpoint_values, aes(x = cutpoint)) +
    geom_line(aes(y = cindex), color = "blue", size = 1) +
    geom_point(aes(y = cindex), color = "blue") +
    
    geom_line(aes(y = neg_log10_p / max(neg_log10_p)),
              color = "red", size = 1) +
    geom_point(aes(y = neg_log10_p / max(neg_log10_p)),
               color = "red") +
    
    scale_y_continuous(
      name = "C-index",
      sec.axis = sec_axis(
        ~ . * max(cutpoint_values$neg_log10_p),
        name = "-log10(p-value)"
      )
    ) +
    
    geom_vline(
      xintercept = opt_cutpoint_output$optimal_cutpoint$cutpoint,
      linetype = "dashed",
      color = "black",
      size = 1
    ) +
    
    theme_minimal() +
    ggtitle("C-index and p-value vs Cutpoint") +
    theme(
      axis.title.y.left = element_text(color = "blue"),
      axis.title.y.right = element_text(color = "red")
    )
  
  return(plot)
}

# Format Cox model output
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
cols <- grep(pattern = "\\.(01|03)[A-Z]$", cols, value = TRUE)
expr.data <- expr.data[, cols]

# Edit column names so they match survival data
colnames(expr.data) <-
  str_replace_all(
    str_remove(colnames(expr.data),
               pattern = "\\.(01|03)[A-Z]$"),
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

# Annotate the genes, leave the ones without HUGO symbol as ENSG
expr.data <- dge.logcpm %>%
  as.data.frame() %>%
  rownames_to_column("id") %>% # Add a new column named "id" with the rownames
  left_join(genes[, c("id", "gene")], by = "id") %>%
  mutate(
    Hugo_Symbol = ifelse(is.na(gene) | gene == "", id, gene),
    Entrez_Gene_Id = id
  ) %>%
  select(-gene, -id)

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
saveRDS(expr.data, file = "survival_results/expr.data_tcga.rds")




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

saveRDS(gsva, "gsva_final.rds")


## RUNNING COX MODELS FOR THERAPEUTIC CLUSTER MARKERS
# Create a prognostic model of survival based on the functional gene expression
# signatures:
expr.data <- readRDS("survival_results/expr.data_tcga.rds")
gsva <- readRDS("survival_results/gsva_alltypes.rds")
gsva <- as.data.frame(t(gsva))
gsva["bcr_patient_barcode"] <- rownames(gsva)

# In case of cluster 9 we only have one gene, so we add the normalized expression 
# of this gene to the gsva table
cluster9 <- expr.data["ARHGDIB", ]
gsva$Cluster09_UP <- cluster9
gsva <- gsva %>%
  relocate(Cluster09_UP, .after = Cluster08_UP)
colnames(gsva)[1:10] <- paste0("Cluster_", 1:10)

# We join the table of unscaled gsva scores to the metadata information.
gsva_metadata <-
  left_join(gsva, metadata, by = "bcr_patient_barcode")

# Adapt the clinical stage
gsva_metadata$tumor_stage <- case_when(
  gsva_metadata$ajcc_pathologic_tumor_stage %in% 
    c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
  gsva_metadata$ajcc_pathologic_tumor_stage %in% 
    c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "Stage II",
  gsva_metadata$ajcc_pathologic_tumor_stage %in% 
    c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Stage III",
  gsva_metadata$ajcc_pathologic_tumor_stage %in% 
    c("Stage IV", "Stage IVA", "Stage IVB") ~ "Stage IV",
  TRUE ~ "Other"
)

gsva_metadata$grade <- case_when(
  gsva_metadata$histological_grade %in% 
    c("G1", "G2", "Low Grade") ~ "Low Grade",
  gsva_metadata$histological_grade %in% 
    c("G3", "G4", "High Grade", "GB") ~ "High Grade",
  gsva_metadata$histological_grade %in%
    c("GX", "[Not Available]", "[Unknown]", "[Discrepancy]") ~ "Unknown",
  TRUE ~ "Other"
)

gsva_metadata$age <- gsva_metadata$age_at_initial_pathologic_diagnosis

# Compute Cox Proportional-Hazards Model per gene set (adjusting by cancer type,
#gender, age and tumor stage).
fix_cov <-
  c("type",
    "gender",
    "age",
    "tumor_stage", 
    "grade")

# Set reference levels
gsva_metadata$type <- relevel(factor(gsva_metadata$type), ref = "THCA")
gsva_metadata$gender <- relevel(factor(gsva_metadata$gender), ref = "MALE")
gsva_metadata$tumor_stage <- relevel(factor(gsva_metadata$tumor_stage), ref = "Stage I")
gsva_metadata$grade <- relevel(factor(gsva_metadata$grade), ref = "Low Grade")


# Cox regression analysis with continuous variable (GSVA score per signature)
# using PFI (Progression Free Survival) as response variable.
signatures <- colnames(gsva_metadata[1:10])
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
saveRDS(models_PFI, "survival_results/cox_models_pfi.rds")
write.table(
  results_PFI,
  file = "survival_results/cox_pfi.tsv",
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
saveRDS(models_OS, "survival_results/cox_models_os.rds")
write.table(
  results_OS,
  file = "survival_results/cox_os.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)


# Plot results from individual Cox models from each cluster into a single plot
forest_plot <- forest_model(model_list = models_PFI, 
                            covariates = signatures,
                            format_options = forest_model_format_options(text_size = 5),
                            merge_models = TRUE,
                            theme = theme_forest() +
                              theme(axis.text.x = element_text(size = 14, color = "black")))

ggsave("survival_results/forest_pfi.png", 
       plot = forest_plot, 
       width = 9, 
       height = 5, 
       dpi = 700)


## Running Cox models independently for each cancer type
# Mapping: TCGA code -> optimal survival endpoint based on cancer type
tcga_survival <- list(
  ACC  = "OS",
  BLCA = "OS",
  BRCA = "PFI",
  CESC = "OS",
  CHOL = "OS",
  COAD = "OS",
  DLBC = "PFI",
  ESCA = "OS",
  GBM  = "OS",
  HNSC = "OS",
  KICH = "OS",
  KIRC = "OS",
  KIRP = "OS",
  LAML = "OS",
  LGG  = "PFI",
  LIHC = "OS",
  LUAD = "OS",
  LUSC = "OS",
  MESO = "OS",
  OV   = "OS",
  PAAD = "OS",
  PCPG = "PFI",
  PRAD = "PFI",
  READ = "PFI",
  SARC = "OS",
  SKCM = "OS",
  STAD = "OS",
  TGCT = "PFI",
  THCA = "PFI",
  THYM = "PFI",
  UCEC = "OS",
  UCS  = "OS",
  UVM  = "OS"
)

# Filter covariates without at least to different values within the specific cancer type
filter_covariates <- function(data, covariates) {
  covariates[sapply(covariates, function(cov) {
    vals <- data[[cov]]
    if (is.factor(vals) || is.character(vals)) {
      length(unique(na.omit(vals))) > 1
    } else {
      TRUE
    }
  })]
}

surv_function <- function(signatures, cancer_type, surv_param){
  gsva_metadata_subset <- gsva_metadata %>%
    filter(type == cancer_type)
  fix_cov <- c("gender",
               "age",
               "tumor_stage",
               "grade")
  fix_cov <- filter_covariates(gsva_metadata_subset, fix_cov)
  formulas <- sapply(signatures,
                     function(x)
                       as.formula(paste0('Surv(', surv_param, ".time, ", surv_param, ')~',
                             paste(c(fix_cov, x), collapse = " + ")
                           )))
  models <-
    lapply(formulas, function(x) {
      coxph(x, data = gsva_metadata_subset)
    })
  
  results <-
    as.data.frame(t(mapply(getResCox_output, models, signatures)))
  return(results)
}

all_results <- list()

for (signature in signatures) {
  signature_results <- list()
  
  for (cancer in unique(gsva_metadata$type)) {
    result <- surv_function(signature, cancer, tcga_survival[cancer])
    
    if (!is.null(result)) {
      # Add cancer type as a column for tracking
      result$cancer_type <- cancer
      result$signature <- signature
      signature_results[[cancer]] <- result
    }
  }
  
  # Combine all cancer-specific results into one data frame
  if (length(signature_results) > 0) {
    combined_results <- bind_rows(signature_results)
    all_results[[signature]] <- combined_results
  }
}

# Save all results to an Excel workbook
write.xlsx(all_results, file = "survival_results/survival_pfi_os_cancer_wise.xlsx")



# Survival curves of cluster 10 signature in each cancer type
for (cancer in unique(gsva_metadata$type)){
  
  # Subset data for this cancer type
  cancer_data <- gsva_metadata %>% filter(type == cancer)
  
  # Get survival time and event column names
  surv_event <- tcga_survival[[cancer]]
  surv_time <- paste0(surv_event, ".time")
  
  # Select optimal cutpoint
  x <- opt_cutpoint(cancer_data, surv_event, cancer, "Cluster_10")
  plot <- opt_cutpoint_plot(x)
  ggsave(paste0("survival_results/KM_curves/cluster10/", cancer, "opt_cutpoint.png"),
         plot = plot,
         width = 10,
         height = 6,
         units = "in",
         dpi = 300)
  cutpoint <- x$optimal_cutpoint$cutpoint
  
  # Assign "Low" and "High" groups based on cutpoint
  cancer_data$Cluster_class <- ifelse(cancer_data$Cluster_10 > cutpoint, "High", "Low")
  cancer_data$Cluster_class <- factor(cancer_data$Cluster_class, levels = c("Low", "High"))
  
  # Create survival object dynamically
  surv_obj <- Surv(time = cancer_data[[surv_time]], event = cancer_data[[surv_event]])
  
  # Get counts per group
  group_counts <- table(cancer_data$Cluster_class)
  legend_labels <- paste0(names(group_counts), " score (n = ", group_counts, ")")
  legend_labels <- legend_labels[match(c("Low", "High"), names(group_counts))]
  
  # Fit and plot KM curve
  fit <- survfit(surv_obj ~ Cluster_class, data = cancer_data)
  survplot <- ggsurvplot(fit,
                         conf.int = TRUE, 
                         risk.table = TRUE, 
                         pval = TRUE,
                         palette = c("#28cfb8", "#FF6663"),
                         legend.labs = legend_labels,
                         title = paste0("Survival - ", cancer, "(", surv_event, ")"))
  print(cancer)
  png(file = paste0("survival_results/KM_curves/cluster10/", cancer, ".png"),
      width = 10,
      height = 6,
      units = "in",
      res = 300)
  print(survplot)
  dev.off()
}





#### RUNNING COX MODELS FOR CLUSTER 4 ENRICHED AMPLIFICATIONS
# Create a prognostic model of survival based on the amplification gene expression
# signatures:
gmt <- getGmt("survival_results/tc4_cnvs.gmt")
genes$id <- gsub("\\.\\d+$", "", genes$id)
ensembl_to_gene <- setNames(genes$gene, genes$id)

# Change from ENSG format to gene names
replace_ids <- function(ids, mapping) {
  vapply(ids, function(id) {
    if (id %in% names(mapping)) mapping[[id]] else id
  }, character(1))
}

updated_sets <- lapply(gmt, function(gene_set) {
  genes <- geneIds(gene_set)
  new_genes <- unique(replace_ids(genes, ensembl_to_gene))
  
  if (length(new_genes) == 0 || all(is.na(new_genes))) {
    message("No mapped genes for: ", 
            setName(gene_set), 
            " — keeping original gene set.")
    return(gene_set)  # Keep the original gene set
  }
  
  GeneSet(
    geneIds = new_genes,
    geneIdType = SymbolIdentifier(),
    setName = setName(gene_set)
  )
})
gsc <- GeneSetCollection(updated_sets)
toGmt(gsc, con = "survival_results/tc4_cnvs_translated.gmt")




gsva <- readRDS("survival_results/gsva_c4_alltypes_translated.rds")
gsva <- as.data.frame(t(gsva))
colnames(gsva) <- gsub("\\.", "-", colnames(gsva))
#colnames(gsva) <- paste0("`", colnames(gsva), "`")
signatures <- colnames(gsva)
gsva["bcr_patient_barcode"] <- rownames(gsva)

# We join the table of unscaled gsva scores to the metadata information.
gsva_metadata <-
  left_join(gsva, metadata, by = "bcr_patient_barcode")

# Adapt the clinical stage
gsva_metadata$tumor_stage <- case_when(
  gsva_metadata$ajcc_pathologic_tumor_stage %in% 
    c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
  gsva_metadata$ajcc_pathologic_tumor_stage %in% 
    c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "Stage II",
  gsva_metadata$ajcc_pathologic_tumor_stage %in% 
    c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Stage III",
  gsva_metadata$ajcc_pathologic_tumor_stage %in% 
    c("Stage IV", "Stage IVA", "Stage IVB") ~ "Stage IV",
  TRUE ~ "Other"
)

gsva_metadata$grade <- case_when(
  gsva_metadata$histological_grade %in% 
    c("G1", "G2", "Low Grade") ~ "Low Grade",
  gsva_metadata$histological_grade %in% 
    c("G3", "G4", "High Grade", "GB") ~ "High Grade",
  gsva_metadata$histological_grade %in%
    c("GX", "[Not Available]", "[Unknown]", "[Discrepancy]") ~ "Unknown",
  TRUE ~ "Other"
)

gsva_metadata$age <- gsva_metadata$age_at_initial_pathologic_diagnosis

# Compute Cox Proportional-Hazards Model per gene set (adjusting by cancer type,
#gender, age and tumor stage).
fix_cov <-
  c("type",
    "gender",
    "age",
    "tumor_stage", 
    "grade")

# Set reference levels
gsva_metadata$type <- relevel(factor(gsva_metadata$type), ref = "THCA")
gsva_metadata$gender <- relevel(factor(gsva_metadata$gender), ref = "MALE")
gsva_metadata$tumor_stage <- relevel(factor(gsva_metadata$tumor_stage), ref = "Stage I")
gsva_metadata$grade <- relevel(factor(gsva_metadata$grade), ref = "Low Grade")


# Cox regression analysis with continuous variable (GSVA score per signature)
# using PFI (Progression Free Survival) as response variable.
formulas_PFI <- sapply(signatures, function(x) {
  varname <- if (!make.names(x) == x)
    paste0("`", x, "`")
  else
    x
  as.formula(paste('Surv(PFI.time, PFI) ~', paste(c(fix_cov, varname), collapse = " + ")))
})

models_PFI <-
  lapply(formulas_PFI, function(x) {
    coxph(x, data = gsva_metadata)
  })

results_PFI <-
  as.data.frame(t(mapply(getResCox_output, models_PFI, paste0("`", signatures, "`"))))
saveRDS(models_PFI, "survival_results/genomic_bands/cox_pfi_genomicbands4.rds")
write.table(
  results_PFI,
  file = "survival_results/genomic_bands/cox_pfi_genomicbands4.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)


# Cox regression analysis with continuous variable (GSVA score per signature)
# using OS (Overall Survival) as response variable.
formulas_OS <- sapply(signatures, function(x) {
  varname <- if (!make.names(x) == x)
    paste0("`", x, "`")
  else
    x
  as.formula(paste('Surv(OS.time, OS) ~', paste(c(fix_cov, varname), collapse = " + ")))
})

models_OS <-
  lapply(formulas_OS, function(x) {
    coxph(x, data = gsva_metadata)
  })

results_OS <-
  as.data.frame(t(mapply(getResCox_output, models_OS, paste0("`", signatures, "`"))))
saveRDS(models_OS, "survival_results/genomic_bands/cox_os_genomicbands4.rds")
write.table(
  results_OS,
  file = "survival_results/genomic_bands/cox_os_genomicbands4.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)


# Plot results from individual Cox models from each cluster into a single plot
forest_plot <- forest_model(model_list = models_OS, 
                            covariates = paste0("`", signatures, "`"),
                            format_options = forest_model_format_options(text_size = 5),
                            merge_models = TRUE,
                            theme = theme_forest() +
                              theme(axis.text.x = element_text(size = 14, color = "black")))

ggsave("survival_results/genomic_bands/forest_os.png", 
       plot = forest_plot, 
       width = 9, 
       height = 5, 
       dpi = 700)

# Cox models of genomic bands per cancer type
all_results <- list()

for (signature in paste0("`", signatures, "`")) {
  signature_results <- list()
  
  for (cancer in unique(gsva_metadata$type)) {
    result <- surv_function(signature, cancer, tcga_survival[cancer])
    
    if (!is.null(result)) {
      # Add cancer type as a column for tracking
      result$cancer_type <- cancer
      result$signature <- signature
      signature_results[[cancer]] <- result
    }
  }
  
  # Combine all cancer-specific results into one data frame
  if (length(signature_results) > 0) {
    combined_results <- bind_rows(signature_results)
    all_results[[signature]] <- combined_results
  }
}

# Save all results to an Excel workbook
write.xlsx(all_results, 
           file = "survival_results/genomic_bands/survival_pfi_os_cancer_wise.xlsx")







# Explore the association of RAP2B expression with survival in ESCA
rap2b <- as.data.frame(t(expr.data["RAP2B",, drop = FALSE]))
rap2b$bcr_patient_barcode <- rownames(rap2b)
rap2b_metadata <-
  left_join(rap2b, gsva_metadata[, grep("^Cluster", colnames(gsva_metadata), value = TRUE, invert = TRUE)], 
            by = "bcr_patient_barcode")
rap2b_metadata_esca <- rap2b_metadata %>%
  filter(type == "ESCA")
res.cox <- coxph(Surv(PFI.time, PFI) ~ type + gender + age + tumor_stage + RAP2B,
                 data = rap2b_metadata)

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
saveRDS(models_PFI, "survival_results/cox_models_pfi.rds")
write.table(
  results_PFI,
  file = "survival_results/cox_pfi.tsv",
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
saveRDS(models_OS, "survival_results/cox_models_os.rds")
write.table(
  results_OS,
  file = "survival_results/cox_os.tsv",
  sep = "\t",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE
)













vars <- gsva[, c("age", "Cluster_10", "cluster2", "cluster3")]


fit <- coxph(Surv(PFI.time, PFI)~ age + Cluster_10, data = x)
vif(fit)


res.cut <- surv_cutpoint(gsva_metadata,
                         time = "PFI.time",
                         event = "PFI",
                         variables = "Cluster_10",
                         minprop = 0.2)

# Calculate 25th and 75th percentiles (Q1 and Q3)
q25 <- quantile(gsva_metadata$Cluster_10, 0.25, na.rm = TRUE)
q75 <- quantile(gsva_metadata$Cluster_10, 0.75, na.rm = TRUE)

# Assign "Low" and "High" groups based on quartiles
gsva_metadata$Cluster10_class[gsva_metadata$Cluster_10 < q25] <- "Low"
gsva_metadata$Cluster10_class[gsva_metadata$Cluster_10 > q75] <- "High"

# Convert to factor (optional, useful for modeling/plotting)
gsva_metadata$Cluster10_class <- factor(gsva_metadata$Cluster10_class, levels = c("Low", "High"))
gsva_metadata$Cluster10_class <- relevel(gsva_metadata$Cluster10_class, ref = "Low")

survplot <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ Cluster10_class, data = gsva_metadata),
                       conf.int = T, risk.table = TRUE, pval = TRUE,
                       palette = c("#28cfb8", "#FF6663"),
                       legend.labs = c("Low score", "High score"))
png(file = "survival_clust10.png",
    width = 10,
    height = 6,
    units = "in",
    res = 300)
survplot
dev.off()





# Survival curve for BRCA
gsva_brca <- gsva_metadata 

res.cox <- coxph(Surv(PFI.time, PFI) ~ type + gender + age + tumor_stage + Cluster10_UP,
                 data = gsva_metadata)

cutpoints <- quantile(gsva_metadata$Cluster10_UP, probs = seq(0.1, 0.9, by = 0.01))

# For each cutpoint, compute C-index and p-value
data <- gsva_filtered[complete.cases(gsva_filtered[, c("time", "event", "type", "gender", "age", "tumor_stage", "Cluster10_UP")]), ]
results <- lapply(cutpoints, function(cut) {
  data$group <- ifelse(data$Cluster10_UP > cut, "high", "low")
  data$time <- data$PFI.time
  data$event <- data$PFI
  surv_obj <- Surv(time = data$time, event = data$event)
  cox_fit <- coxph(Surv(time = data$time, event = data$event) ~ type + gender + age + tumor_stage + group,
                   data = data)
  
  # Extract the model frame used internally (with complete cases only)
  pred <- predict(cox_fit)
  used_data <- model.frame(cox_fit)
  
  
  # Compute C-index
  library(survcomp)
  c_index <- concordance.index(pred, surv.time = data$PFI.time, 
                               surv.event = data$PFI)$c.index
  
  # Extract p-value
  pval <- summary(cox_fit)$coefficients["grouplow", "Pr(>|z|)"]
  
  return(data.frame(cut = cut, c_index = c_index, pval = pval))
})

results_df <- do.call(rbind, results)

# Select optimal cutpoint (e.g. max C-index with p < 0.05)
valid_results <- results_df %>% filter(pval < 0.05)
optimal_row <- valid_results[which.max(valid_results$c_index), ]

# Clean p-values: clip to avoid extreme scaling (optional)
results_df$pval <- pmin(results_df$pval, 1)
results_df$pval[results_df$pval < 1e-10] <- 1e-10

# Define desired axis limits
cindex_min <- 0.5
cindex_max <- 1.0
pval_min <- 0.0
pval_max <- 0.05

# Scale factor between axes based on desired limits
scale_factor <- (cindex_max - cindex_min) / (pval_max - pval_min)

# Plot
ggplot(results_df, aes(x = cut)) +
  # C-index curve
  geom_line(aes(y = c_index, color = "C-index"), size = 1) +
  
  # p-value scaled
  geom_line(aes(y = pval * scale_factor + cindex_min, color = "p-value"),
            linetype = "dashed", size = 1) +
  
  # Optimal cutpoint vertical line
  geom_vline(xintercept = optimal_row$cut, linetype = "dashed", color = "darkgreen") +
  
  # Axis definitions
  scale_y_continuous(
    name = "C-index",
    limits = c(cindex_min, cindex_max),
    sec.axis = sec_axis(
      trans = ~ (. - cindex_min) / scale_factor,
      name = "p-value",
      breaks = seq(0, 1, 0.2)
    )
  ) +
  
  # Colors and theme
  scale_color_manual(values = c("C-index" = "blue", "p-value" = "red")) +
  labs(
    title = "Optimal Cutpoint Selection",
    x = "Cutpoint",
    color = "Metric"
  ) +
  theme_minimal() +
  theme(
    axis.title.y.left = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "red"),
    legend.position = "top"
  )

gsva_brca <- gsva_metadata %>%
  filter(type == "BRCA")
res.cut <- surv_cutpoint(gsva_metadata,
                         time = "PFI.time",
                         event = "PFI",
                         variables = "Cluster10_UP",
                         minprop = 0.2)

# Calculate 25th and 75th percentiles (Q1 and Q3)
q25 <- quantile(gsva_brca$Cluster10_top50_UP, 0.25, na.rm = TRUE)
q75 <- quantile(gsva_brca$Cluster10_top50_UP, 0.75, na.rm = TRUE)

# Assign "Low" and "High" groups based on quartiles
gsva_metadata$Cluster10_class[gsva_metadata$Cluster10_UP < res.cut$cutpoint$cutpoint] <- "Low"
gsva_metadata$Cluster10_class[gsva_metadata$Cluster10_UP > res.cut$cutpoint$cutpoint] <- "High"

# Convert to factor (optional, useful for modeling/plotting)
gsva_metadata$Cluster10_class <- factor(gsva_metadata$Cluster10_class, levels = c("Low", "High"))
gsva_metadata$Cluster10_class <- relevel(gsva_metadata$Cluster10_class, ref = "Low")

survplot <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ Cluster10_class, data = gsva_metadata),
                 conf.int = T, risk.table = TRUE, pval = TRUE,
                 palette = c("#28cfb8", "#FF6663"),
                 legend.labs = c("Low score", "High score"))
png(file = "survival_clust10.png",
    width = 6,
    height = 6,
    units = "in",
    res = 300)
survplot
dev.off()


# BRCA
gsva_brca <- gsva_metadata %>%
  filter(type == "COAD")
res.cut <- surv_cutpoint(gsva_brca,
                         time = "PFI.time",
                         event = "PFI",
                         variables = "Cluster10_UP",
                         minprop = 0.2)

# Calculate 25th and 75th percentiles (Q1 and Q3)
q25 <- quantile(gsva_brca$Cluster10_top50_UP, 0.25, na.rm = TRUE)
q75 <- quantile(gsva_brca$Cluster10_top50_UP, 0.75, na.rm = TRUE)

# Assign "Low" and "High" groups based on quartiles
gsva_brca$Cluster10_class[gsva_brca$Cluster10_UP < res.cut$cutpoint$cutpoint] <- "Low"
gsva_brca$Cluster10_class[gsva_brca$Cluster10_UP > res.cut$cutpoint$cutpoint] <- "High"

# Convert to factor (optional, useful for modeling/plotting)
gsva_brca$Cluster10_class <- factor(gsva_brca$Cluster10_class, levels = c("Low", "High"))
gsva_brca$Cluster10_class <- relevel(gsva_brca$Cluster10_class, ref = "Low")

survplot <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ Cluster10_class, data = gsva_brca),
                       conf.int = T, risk.table = TRUE, pval = TRUE,
                       palette = c("#28cfb8", "#FF6663"),
                       legend.labs = c("Low score", "High score"))
png(file = "survival_clust10_BRCA.png",
    width = 8,
    height = 6,
    units = "in",
    res = 300)
survplot
dev.off()


gsva <- gsva(
  expr = expr.data,
  gset.idx.list = list(c4_specific = c("KRT17", "KRT16", "KRT14", "JUP", "KRT19")),
  method = "ssgsea",
  ssgsea.norm = TRUE,
  min.sz = 5,
  max.sz = 1800,
  parallel.sz = 5
)









# Survival curves for Cluster 10 signature for PAAD and HSC
subset_data <- gsva_metadata
median_val <- median(gsva_metadata$c4_specific, na.rm = TRUE)

# Classify based on median
gsva_metadata$clust10_pfi <- ifelse(gsva_metadata$c4_specific < median_val, "Low", "High")
gsva_metadata$clust10_pfi <- factor(gsva_metadata$clust10_pfi)


subset_data <- gsva_metadata

# Calculate 25th and 75th percentiles (Q1 and Q3)
q25 <- quantile(subset_data$c4_specific, 0.25, na.rm = TRUE)
q75 <- quantile(subset_data$c4_specific, 0.75, na.rm = TRUE)


# Assign "Low" and "High" groups based on quartiles
subset_data$clust4_os[subset_data$c4_specific < q25] <- "Low"
subset_data$clust4_os[subset_data$c4_specific > q75] <- "High"

# Convert to factor (optional, useful for modeling/plotting)
subset_data$clust4_os <- factor(subset_data$clust4_os, levels = c("Low", "High"))

# Set the reference levels for each categorical variable. We try to use as reference 
# levels the ones with lower malignancy.
subset_data$clust4_os <- relevel(subset_data$clust4_os, ref = "Low")

ggsurv1 <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ clust4_os, data = subset_data),
                      conf.int = T, risk.table = TRUE,
                      palette = c("#28cfb8", "#FF6663"),
                      legend.labs = c("Low score", "High score"))

png("clust10_pfi_curve.png", width = 8, height = 7, units = "in", res = 500)
ggsurv1
dev.off()



res.cut <- surv_cutpoint(gsva_metadata,
                         time = "OS.time", event = "OS",
                         variables = "Cluster04_UP")
gsva_metadata$clust10_pfi <- ifelse(gsva_metadata$Cluster04_UP < res.cut$cutpoint$cutpoint, "Low", "High")
gsva_metadata$clust10_pfi <- factor(gsva_metadata$clust10_pfi)

# Set the reference levels for each categorical variable. We try to use as reference 
# levels the ones with lower malignancy.
gsva_metadata$clust10_pfi<- relevel(gsva_metadata$clust10_pfi, ref = "Low")

ggsurv1 <- ggsurvplot(survfit(Surv(PFI.time, PFI) ~ clust10_pfi, data = gsva_metadata),
                      conf.int = T, risk.table = TRUE,
                      palette = c("#28cfb8", "#FF6663"),
                      legend.labs = c("Low score", "High score"))

png("clust04_os_curve.png", width = 8, height = 7, units = "in", res = 500)
ggsurv1
dev.off()



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
gene_set <- collection@.Data[[10]]@geneIds
gene_sets <- lapply(gene_set, function(gene) setdiff(gene_set, gene))
names(gene_sets) <- paste("cluster10", gene_set, sep = "_")
names(gene_sets) <- gsub("-", "_", names(gene_sets))

# Also include the full signature if you want to compare
gene_sets[["original_signature"]] <- gene_set

# Run GSVA (or ssGSEA) on all signatures in one shot
gsva <- gsva(
  expr = expr.data,
  gset.idx.list = gene_sets,
  method = "ssgsea",
  ssgsea.norm = TRUE,
  min.sz = 3,
  max.sz = 1800,
  parallel.sz = 5
)

gsva <- as.data.frame(t(gsva))
signatures <- colnames(gsva)
gsva["bcr_patient_barcode"] <- rownames(gsva)
gsva_metadata <-
  left_join(gsva, metadata, by = "bcr_patient_barcode")

library(survival)

results <- data.frame(
  signature = character(),
  HR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  pvalue = numeric(),
  stringsAsFactors = FALSE
)

fix_cov <-
  c("type",
    "gender",
    "age_at_initial_pathologic_diagnosis",
    "clinical_stage")

for (sig in signatures) {
  formula <- as.formula(paste(
                        'Surv(PFI.time, PFI)~',
                        paste(c(fix_cov, sig), collapse = " + ")))
  model <- coxph(formula,
                 data = gsva_metadata)
  s <- summary(model)
  
  results <- rbind(results, data.frame(
    signature = sig,
    HR = s$coefficients[sig, "exp(coef)"],
    CI_lower = s$conf.int[sig, "lower .95"],
    CI_upper = s$conf.int[sig, "upper .95"],
    pvalue = s$coefficients[sig, "Pr(>|z|)"]
  ))
}
