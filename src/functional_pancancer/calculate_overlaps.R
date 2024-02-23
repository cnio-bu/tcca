library(fgsea)
library(tidyverse)

gsets <- gmtPathways("reference/combined_gsets_functional.gmt")

gsets <- lapply(gsets, FUN = function(x){x[x != ""]})

jaccards <- unlist(lapply(combn(gsets, 2, simplify = FALSE), function(x) {
    length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))

gsets_comps <- combn(names(gsets), 2, simplify = FALSE)
gsets_comps <- as.data.frame(gsets_comps)
gsets_comps <- as.data.frame(t(gsets_comps))
gsets_comps$comp <- paste0(gsets_comps$V1, "_VS_", gsets_comps$V2)


jaccards <- as.data.frame(jaccards)
jaccards$comparison <- gsets_comps$comp

jaccards <- jaccards %>%
    separate(comparison, into = c("gset1", "gset2"), sep = "_VS_") %>%
    mutate(
        jaccards = round(jaccards, digits = 3)
    )
    

write.table(jaccards, file = "results/functional/jaccard_mat_long.tsv")

