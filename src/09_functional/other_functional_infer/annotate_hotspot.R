library(dplyr)
library(tidyverse)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(fgsea)
library(openxlsx)
setwd("/home/lmgonzalezb/Documents/bc-meta/functional_mps/hotspot/")
gene_modules <- read.table("gene_modules.tsv", sep = "\t", header = TRUE)

# Load gene set collection from MSigDB
H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
C2 <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  filter(gs_subcat != "CGP") %>% 
  dplyr::select(gs_name, gene_symbol)
C5_BP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = c("GO:BP")) %>% 
  dplyr::select(gs_name, gene_symbol)
C5_CC <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = c("GO:CC")) %>% 
  dplyr::select(gs_name, gene_symbol)
C5_MF <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = c("GO:MF")) %>% 
  dplyr::select(gs_name, gene_symbol)
C6 <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, gene_symbol)
gene_sets <- rbind(H, C2, C5_BP, C5_CC, C5_MF, C6)

pathway.list <- gene_sets %>%
  group_by(gs_name) %>%
  summarise(Genes = list(gene_symbol), .groups = 'drop') %>%
  deframe()

# Create gene modules list
gene_modules.list <- list()
modules <- sort(na.omit(unique(gene_modules$Module)))
for (num in modules){
  genes <- gene_modules %>% 
    filter(Module == num) %>% 
    pull(Gene)
  gene_modules.list[[paste0("Module_", num)]] <- genes
}

# Compute overrepresentation of gene modules across MSigDB collections
func_annot <- lapply(gene_modules.list, function(program) {
  fgRes <- fgsea::fora(pathways = pathway.list,
                       genes = program,
                       universe = as.vector(gene_modules$Gene))
  
  fgRes <- fgRes[fgRes$pval <= 0.05,]
  return(fgRes)
})

# Save results as a excel file
wb <- createWorkbook()
# Write each data frame to a separate sheet
for (i in seq_along(func_annot)) {
  addWorksheet(wb, sheetName = names(func_annot)[i]) # Use the list names as sheet names
  writeData(wb, sheet = names(func_annot)[i], x = func_annot[[i]])
}

# Save the workbook to a file
saveWorkbook(wb, file = "fgsea_msigdb_ora.xlsx", overwrite = TRUE)



# Compute overrepresentation of gene modules across our custom functional gene sets
gmt.df <- read.gmt("../../bc-meta_repo/bc-meta/reference/combined_gsets_functional.gmt")
gmt.list <- split(x = gmt.df$gene, f = gmt.df$term)
gmt.list <- lapply(gmt.list, function(x) x[x != ""])

func_annot <- lapply(gene_modules.list, function(program) {
  fgRes <- fgsea::fora(pathways = gmt.list,
                       genes = program,
                       universe = as.vector(gene_modules$Gene))
  
  fgRes <- fgRes[fgRes$pval <= 0.05,]
  return(fgRes)
})

wb <- createWorkbook()
for (i in seq_along(func_annot)) {
  addWorksheet(wb, sheetName = names(func_annot)[i]) # Use the list names as sheet names
  writeData(wb, sheet = names(func_annot)[i], x = func_annot[[i]])
}
saveWorkbook(wb, file = "fgsea_custom_genesets_ora.xlsx", overwrite = TRUE)

