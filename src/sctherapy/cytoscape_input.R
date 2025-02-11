library(dplyr)
library(tidyverse)

### FUNCTIONS ###
jaccard_similarity <- function(vector1, vector2) {
  intersection <- length(intersect(vector1, vector2))
  union <- length(union(vector1, vector2))
  return(intersection / union)
}

setwd("/home/lmgonzalezb/Documents/bc-meta/sctherapy")

data <- read.table("full_table_drug_prediction.tsv")
data$study_sample <- paste0(sub("\\..*", "", data$Subclone), "_", data$Sample)

drug_subclone <- data %>%
  select(Subclone, Drug_Name) %>%
  distinct()

# Create a list with drug predictions per subclone
subclones_drug_list <- lapply(subclones, function(subclone) {
  drugs <- drug_subclone$Drug_Name[drug_subclone$Subclone == subclone]
  return(drugs)
})
subclones <- unique(data$Subclone)
names(subclones_drug_list) <- subclones


# Compute shared drugs between subclones
pairs <- t(combn(subclones, 2))
pairs <- as.data.frame(pairs)
colnames(pairs) <- c("Source", "Target")

edge_table <- pairs %>%
  rowwise() %>%
  mutate(Weight = jaccard_similarity(subclones_drug_list[[Source]], 
                                     subclones_drug_list[[Target]]))

write.table(edge_table, "network_edges.tsv", row.names = FALSE)
