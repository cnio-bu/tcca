library(beyondcell)
library(dplyr)

setwd("/Users/mariagb/Documents/bc_meta/beyondcell")

ssc <- GetCollection(SSc, include.pathways = FALSE)

# Function to convert GeneSet object from Beyondcell to GMT format
geneset_to_gmt <- function(geneset_obj, output_file = "signatures.gmt", 
                           include_drug_info = TRUE) {
  
  # Extract list of gene sets
  genelist <- geneset_obj@genelist
  
  # Extract drug information if available
  has_info <- nrow(geneset_obj@info) > 0
  
  if (has_info && include_drug_info) {
    # Group information by IDs to obtain unique descriptions
    info_df <- geneset_obj@info
    sig_info <- info_df %>%
      dplyr::group_by(IDs) %>%
      dplyr::summarise(
        description = paste0(
          unique(drugs), " | ",
          unique(MoAs), " | ",
          unique(targets), " | ",
          unique(studies)
        ) %>% paste(collapse = "; "),
        .groups = "drop"
      )
  }
  
  # Open file connection
  con <- file(output_file, "w")
  
  genesets_written <- 0
  
  # Iterate over each gene set
  for (geneset_name in names(genelist)) {
    
    # Get up and down genes
    genes_up <- genelist[[geneset_name]]$up
    genes_down <- genelist[[geneset_name]]$down
    
    # Get description if exists
    if (has_info && include_drug_info) {
      desc_row <- sig_info[sig_info$IDs == geneset_name, ]
      if (nrow(desc_row) > 0) {
        description <- desc_row$description[1]
      } else {
        description <- "ssc"
      }
    } else {
      description <- "ssc"
    }
    
    # Write down UP geneset if exists
    if (!is.null(genes_up) && length(genes_up) > 0) {
      geneset_name_up <- paste0(geneset_name, "_UP")
      line <- paste(c(geneset_name_up, description, genes_up), collapse = "\t")
      writeLines(line, con)
      genesets_written <- genesets_written + 1
    }
    
    # Write down DOWN geneset if exists
    if (!is.null(genes_down) && length(genes_down) > 0) {
      geneset_name_down <- paste0(geneset_name, "_DOWN")
      line <- paste(c(geneset_name_down, description, genes_down), collapse = "\t")
      writeLines(line, con)
      genesets_written <- genesets_written + 1
    }
  }
  
  # Close connection
  close(con)
  
  cat("GMT file saved to:", output_file, "\n")
  cat("Total genesets written:", genesets_written, "\n")
  if (has_info && include_drug_info) {
    cat("Drug information included in descriptions\n")
  }
}


# Get gmt
geneset_to_gmt(ssc, "SSc.gmt", include_drug_info = FALSE)

# Join the SSc collection gmt file to the immunotherapy signatures
merge_gmt_files <- function(gmt_file1, gmt_file2, output_file) {
  
  lines1 <- readLines(gmt_file1)
  lines2 <- readLines(gmt_file2)
  
  all_lines <- c(lines1, lines2)
  
  writeLines(all_lines, output_file)
  
  cat("Merged GMT file saved:", output_file, "\n")
  cat("File 1 genesets:", length(lines1), "\n")
  cat("File 2 genesets:", length(lines2), "\n")
  cat("Total genesets:", length(all_lines), "\n")
}

# Combinar dos GMT
merge_gmt_files(
  gmt_file1 = "SSc.gmt",
  gmt_file2 = "immunotherapy.gmt",
  output_file = "bc_immuno.gmt"
)
