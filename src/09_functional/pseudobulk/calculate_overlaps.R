library(fgsea)
library(tidyverse)

gsets <- gmtPathways("reference/combined_gsets_functional.gmt")

gsets <- lapply(gsets, FUN = function(x){x[x != ""]})

## Filter gsets with < 15 members
gsets <- gsets[lapply(gsets, length) >= 15]

jaccards <- unlist(lapply(combn(gsets, 2, simplify = FALSE), function(x) {
    length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))

intersections <- unlist(lapply(combn(gsets, 2, simplify = FALSE), function(x) {
    length(intersect(x[[1]], x[[2]]))}))

gsets_comps <- combn(names(gsets), 2, simplify = FALSE)
gsets_comps <- as.data.frame(gsets_comps)
gsets_comps <- as.data.frame(t(gsets_comps))
gsets_comps$comp <- paste0(gsets_comps$V1, "_VS_", gsets_comps$V2)


jaccards <- as.data.frame(jaccards)
jaccards$comparison <- gsets_comps$comp
jaccards$intersection <- intersections

jaccards <- jaccards %>%
    separate(comparison, into = c("gset1", "gset2"), sep = "_VS_") %>%
    mutate(
        jaccards = round(jaccards, digits = 3),
    ) %>%
    rowwise() %>%
    mutate(
        gset1_size = length(gsets[[gset1]]),
        gset2_size = length(gsets[[gset2]])
    )
    
jaccards_prop <- jaccards %>%
    mutate(
        gset1_prop = intersection / gset1_size,
        gset2_prop = intersection / gset2_size,
        to_check = (gset1_prop >= 0.5 | gset2_prop >= 0.5)
    )

jacs_to_remove <- jaccards_prop %>%
    filter(to_check) %>%
    rowwise() %>%
    mutate(
        "gset_to_keep" = case_when(
            gset1_size <= gset2_size ~ "gset1",
            gset2_size < gset1_size ~ "gset2"
         )
    )

write.table(jaccards, file = "results/functional/jaccard_mat_long.tsv")

