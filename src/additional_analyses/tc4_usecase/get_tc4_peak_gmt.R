library(tidyverse)

## FUNCTIONS
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

read_gmt_list <- function(gmt) {
  gmt_lines <- readLines(gmt)
  
  signatures <- lapply(gmt_lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    genes <- parts[-c(1, 2)]  # Remove name and description
    return(genes)
  })
}


cnvs_tc <- read.table("/local/bc_meta/scevan/cytobands_CNV/top_amp_del_therapeutic_clusters.tsv", header = T, sep = "\t") %>%
  filter(cluster == 4, 
         arm == "chr3q" | cnv %in% c("a_chr1q21.3", "a_chr17q21.2"),
         type == "a")

all_cnv_genes <- list()
combined_genes <- c()
cnvs_chr3 <- cnvs_tc$cnv[grepl("^a_chr3", cnvs_tc$cnv)]

for (cnv in  cnvs_tc$cnv){
  cnv_name <- gsub("a_", "", cnv)
  cnv_name <- gsub("chr", "", cnv_name)
  
  cnv_genes <- get_genes_by_cytoband(cnv_name)
  cnv_genes <- cnv_genes %>%
    mutate(gene = ifelse(external_gene_name != "", external_gene_name, ensembl_gene_id)) %>%
    dplyr::select(gene) %>%
    unlist() %>%
    as.character() %>%
    unique()
  
  if (cnv %in% cnvs_chr3) {
    combined_genes <- c(combined_genes, cnv_genes) %>% unique()
  }
  
  all_cnv_genes[[cnv_name]] <- cnv_genes
}

all_cnv_genes[["3q25-q29"]] <- combined_genes

## Now, add signatures of the intersection with top TC4 markers
markers <- read_gmt_list("genes_by_cluster.gmt")
markers_tc4 <- markers[[4]]

for (cnv in names(all_cnv_genes)) {
  intersec <- intersect(all_cnv_genes[[cnv]], markers_tc4)
  all_cnv_genes[[paste0(cnv, "_markers")]] <- intersec
}

## Open writing connection, write genes and close
con <- file("gdsc/tc4_cnvs.gmt", open = "wt")

for (cnv in names(all_cnv_genes)) {
  genes <- all_cnv_genes[[cnv]]
  line <- c(cnv, "na", genes)
  writeLines(paste(line, collapse = "\t"), con)
}

close(con)
