library(data.table)

setwd("/storage/scratch01/shared/projects/bc-meta/")

out.dir <- "functional_nmf/therapeutic_from_mps/mp_high_cells"

frac_min <- 0.40

# Get MP-high cells
files <- list.files(out.dir, pattern = "^mp_.*\\.rds$", full.names = TRUE)
mp_all <- rbindlist(lapply(files, readRDS))

mp_all[, padj := p.adjust(pval, method = "BH")]
mp_all[, is_high := padj < 0.05 & frac_detected >= frac_min]

saveRDS(mp_all, file.path(out.dir, "mp_scores_all.rds"))

summary_mp <- mp_all[, .(
  n_cells = .N,
  n_high = sum(is_high),
  pct_high = round(100 * mean(is_high), 1),
  n_samples_with_high = uniqueN(sample[is_high])
), by = mp][order(-pct_high)]

print(summary_mp)

mp_high_cells <- mp_all[is_high == TRUE]
mp_groups <- mp_high_cells[, .N, by = .(sample, mp)][N >= 20]

saveRDS(mp_high_cells, file.path(out.dir, "mp_high_cells.rds"))
saveRDS(mp_groups, file.path(out.dir, "mp_groups_viable.rds"))

cat("Groups MP × sample:", nrow(mp_groups), "\n")
cat("Represented MPs:", uniqueN(mp_groups$mp), "\n")
cat("Represented samples:", uniqueN(mp_groups$sample), "\n")