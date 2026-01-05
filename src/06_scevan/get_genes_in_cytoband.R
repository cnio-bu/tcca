get_genes_by_cytoband <- function(cytobands, dataset = "hsapiens_gene_ensembl") {
  library(biomaRt)
  library(dplyr)
  library(stringr)
  
  ensembl <- useEnsembl(biomart = "genes", dataset = dataset)
  
  # Obtener posiciones de bandas
  cytoband_table <- getBM(
    attributes = c("chromosome_name", "band", "start_position", "end_position"),
    mart = ensembl
  ) %>%
    mutate(cytoband = paste0(chromosome_name, band)) %>%
    filter(cytoband %in% cytobands)
  
  if (nrow(cytoband_table) == 0) {
    warning("Ninguna de las bandas especificadas fue encontrada.")
    return(data.frame())
  }
  
  # Obtener todos los genes con posiciones
  all_genes <- getBM(
    attributes = c("external_gene_name", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "gene_biotype"),
    mart = ensembl
  )
  
  # Cruzar genes con bandas (más eficiente)
  results <- list()
  for (i in seq_len(nrow(cytoband_table))) {
    band <- cytoband_table[i, ]
    genes_in_band <- all_genes %>%
      filter(
        chromosome_name == band$chromosome_name,
        start_position <= band$end_position,
        end_position >= band$start_position
      ) %>%
      mutate(cytoband = band$cytoband)
    results[[i]] <- genes_in_band
  }
  
  result_table <- bind_rows(results)
  return(result_table)
}




cnv_specific_7m_named <- gsub("a_", "", cnv_specific_7m)
cnv_specific_7m_named <- gsub("d_", "", cnv_specific_7m_named)
cnv_specific_7m_named <- gsub("chr", "", cnv_specific_7m_named)

cnv_specific_7m_genes <- get_genes_by_cytoband(cnv_specific_7m_named)

setwd("/local/renacer/")
genes_en_bandas <- get_genes_by_cytoband(c("17q21.2"))
write.table(genes_en_bandas, "genes_in_17q21.2.tsv", sep = "\t")
