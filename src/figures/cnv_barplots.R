library(GenVisR)
library(tidyverse)
library(ggplot2)
library(GenomicRanges)
library(patchwork)
library(AnnotationHub)

## Load CNV data
setwd("/local/bc_meta/scevan/cytobands_CNV/")
cs <- read.table("/local/bc_meta/scevan/cnv_segments_clones_lvl2_cytobands.tsv", header = T, sep = "\t")

## This is a table of averages of amplification and deletions of each region
amp <- apply(cs, 1, function(row) mean(row[row > 2], na.rm = TRUE))
del <- apply(cs, 1, function(row) mean(row[row < 2], na.rm = TRUE))
amp.del <- cbind(amp, del)
amp.del <- replace(amp.del, is.nan(amp.del), 0)

## Make a granges of the cytobands of the human genome
proxy <- "lserranor@cnio.es"
proxy <- httr::use_proxy(Sys.getenv('http_proxy'))
httr::set_config(proxy)
AnnotationHub::setAnnotationHubOption("PROXY", proxy)
AnnotationHub::getAnnotationHubOption("LOCAL")

hub <- AnnotationHub()

hub_hg38 <- subset(hub, 
                   (hub$species == "Homo sapiens") & (hub$genome == "hg38"))

cytobands  <- hub_hg38[[797]]
cytobands$custom_name <- paste0(seqnames(cytobands), cytobands$name)

dfA <- data.frame(names = cytobands$custom_name,
                  chrom = as.data.frame(seqnames(cytobands)))
dfB <- as.data.frame(ranges(cytobands))

cn.data <- cbind(dfA, dfB)
cn.data <- cn.data[,-5]

## Join cytobands to CNV data and set long format
cs <- rownames_to_column(as.data.frame(cs), var = "names")
cs <- merge(cn.data, cs, by = "names")

cs.long <- cs %>%
  pivot_longer(-c(names, value, start, end), names_to = "sampleID", values_to = "mean")

cs.long <- data.frame(region_name = cs.long$names,
                      sample = cs.long$sampleID,
                      chromosome = cs.long$value,
                      start = cs.long$start,
                      end = cs.long$end,
                      probes = 1,
                      segmean = cs.long$mean)

cs.long <- replace(cs.long, is.na(cs.long), 2) %>%
  mutate(sample = str_replace_all(sample, "\\.", "-")) # This reverts scevan name changing

cs.long$sample <- gsub("GSM5645908_Breast_1_biol-rep", "GSM5645908_Breast_1_biol.rep", cs.long$sample) # This three fix three exceptions
cs.long$sample <- gsub("SyS11-met", "SyS11.met", cs.long$sample)
cs.long$sample <- gsub("Travaglini_Krasnow_2020_distal-1b", "Travaglini_Krasnow_2020_distal 1b", cs.long$sample)


## Calculate ITH 
samples <- strsplit(unique(cs.long$sample), "__")
samples <- as.data.frame(do.call(rbind, samples))
rownames(samples) <- unique(cs.long$sample)
samples <- rownames_to_column(samples, var = "sample")
colnames(samples) <- c("sample","study", "sample_original", "clone")
samples <- mutate(samples,
                  study__sample = paste(study, sample_original, sep = "__"))

ith <- as.data.frame(table(samples$study__sample))
ith <- separate(ith, Var1, into = c("study", "sample"), sep = "__", remove = F)
colnames(ith) <- c("study__sample", "study", "sample", "nclones")
write.table(ith, "/local/bc_meta/scevan/genomic_ith.tsv", row.names = F, sep = "\t")

## Load clinical and TCs metadata
clinical <- read.table("/local/bc_meta/clinical_metadata_v4_clean.tsv", header = T, sep = "\t")
clinical <- clinical %>%
  mutate(true_patient = paste0(study, "_", patient),
         study__sample = paste(study, sample, sep = "__")) %>%
  filter(study__sample %in% samples$study__sample) # Keep samples for which we have CNV data

metacom <- readRDS("/local/bc_meta/scevan/threapeutic_clusters.rds")
metacom <- metacom %>%
  mutate(sample = str_replace_all(subclone, "\\.", "__") %>% str_replace("(\\d+)$", "subclone\\1"))
metacom$sample <- gsub("GSM5645908_Breast_1_biol__rep", "GSM5645908_Breast_1_biol.rep", metacom$sample) # Fix manual exceptions
metacom$sample <- gsub("SyS11__met", "SyS11.met", metacom$sample)

## Add samples IDs and other metadata to cs.long
cs.long.merged <- left_join(cs.long, samples, by = "sample")
cs.long.merged <- left_join(cs.long.merged, metacom, by = "sample")
cs.long.merged <- left_join(cs.long.merged, 
                            dplyr::select(clinical, study__sample, tumor_type, stage, sample_type, treated, refined_tumor_site), 
                            by = "study__sample")


## FIGURE 1: Global (all, no cell lines, only primary, only met) #################
## Create layers for well known oncogenes
oncogenes <- read.table("/local/bc_meta/scevan/oncogenes_coordinates.tsv", header = T, sep = "\t")

layer1 <- geom_vline(data = cs.long.merged %>% filter(chromosome == oncogenes[1,]$chr),
                     aes(xintercept = oncogenes[1,]$middle), col = "#3a3b3c", linetype = 2, size = 1)
layer2 <- geom_vline(data = cs.long.merged %>% filter(chromosome == oncogenes[2,]$chr),
                     aes(xintercept = oncogenes[2,]$middle), col = "#3a3b3c", linetype = 2, size = 1)
layer3 <- geom_vline(data = cs.long.merged %>% filter(chromosome == oncogenes[3,]$chr),
                     aes(xintercept = oncogenes[3,]$middle), col = "#3a3b3c", linetype = 2, size = 1)
layer4 <- geom_vline(data = cs.long.merged %>% filter(chromosome == oncogenes[4,]$chr),
                     aes(xintercept = oncogenes[4,]$middle), col = "#3a3b3c", linetype = 2, size = 1)
layer5 <- geom_vline(data = cs.long.merged %>% filter(chromosome == oncogenes[5,]$chr),
                     aes(xintercept = oncogenes[5,]$middle), col = "#3a3b3c", linetype = 2, size = 1)
layer6 <- geom_vline(data = cs.long.merged %>% filter(chromosome == oncogenes[6,]$chr),
                     aes(xintercept = oncogenes[6,]$middle), col = "#3a3b3c", linetype = 2, size = 1)
layer7 <- geom_vline(data = cs.long.merged %>% filter(chromosome == oncogenes[7,]$chr),
                     aes(xintercept = oncogenes[7,]$middle), col = "#3a3b3c", linetype = 2, size = 1)
layer8 <- geom_vline(data = cs.long.merged %>% filter(chromosome == oncogenes[8,]$chr),
                     aes(xintercept = oncogenes[8,]$middle), col = "#3a3b3c", linetype = 2, size = 1)
layer9 <- geom_vline(data = cs.long.merged %>% filter(chromosome == oncogenes[9,]$chr),
                     aes(xintercept = oncogenes[9,]$middle), col = "#3a3b3c", linetype = 2, size = 1)
layer10 <- geom_vline(data = cs.long.merged %>% filter(chromosome == oncogenes[10,]$chr),
                      aes(xintercept = oncogenes[10,]$middle), col = "#3a3b3c", linetype = 2, size = 1)
layer11 <- geom_vline(data = cs.long.merged %>% filter(chromosome == oncogenes[11,]$chr),
                      aes(xintercept = oncogenes[11,]$middle), col = "#3a3b3c", linetype = 2, size = 1)

layer1.2 <- geom_text(data = cs.long.merged %>% filter(chromosome == oncogenes[1,]$chr),
                      aes(x=(oncogenes[1,]$middle+18500000), y=0.75), label=oncogenes[1,]$Gene, angle = 90, size = 9)
layer2.2 <- geom_text(data = cs.long.merged %>% filter(chromosome == oncogenes[2,]$chr),
                      aes(x=(oncogenes[2,]$middle+18500000), y=0.75), label=oncogenes[2,]$Gene, angle = 90, size = 9)
layer3.2 <- geom_text(data = cs.long.merged %>% filter(chromosome == oncogenes[3,]$chr),
                      aes(x=(oncogenes[3,]$middle+18500000), y=0.75), label=oncogenes[3,]$Gene, angle = 90, size = 9)
layer4.2 <- geom_text(data = cs.long.merged %>% filter(chromosome == oncogenes[4,]$chr),
                      aes(x=(oncogenes[4,]$middle-18500000), y=0.75), label=oncogenes[4,]$Gene, angle = 90, size = 9)
layer5.2 <- geom_text(data = cs.long.merged %>% filter(chromosome == oncogenes[5,]$chr),
                      aes(x=(oncogenes[5,]$middle+18500000), y=0.75), label=oncogenes[5,]$Gene, angle = 90, size = 9)
layer6.2 <- geom_text(data = cs.long.merged %>% filter(chromosome == oncogenes[6,]$chr),
                      aes(x=(oncogenes[6,]$middle+18500000), y=0.75), label=oncogenes[6,]$Gene, angle = 90, size = 9)
layer7.2 <- geom_text(data = cs.long.merged %>% filter(chromosome == oncogenes[7,]$chr),
                      aes(x=(oncogenes[7,]$middle+18500000), y=0.90), label=oncogenes[7,]$Gene, angle = 90, size = 9)
layer8.2 <- geom_text(data = cs.long.merged %>% filter(chromosome == oncogenes[8,]$chr),
                      aes(x=(oncogenes[8,]$middle-18500000), y=0.75), label=oncogenes[8,]$Gene, angle = 90, size = 9)
layer9.2 <- geom_text(data = cs.long.merged %>% filter(chromosome == oncogenes[9,]$chr),
                      aes(x=(oncogenes[9,]$middle+18500000), y=0.75), label=oncogenes[9,]$Gene, angle = 90, size = 9)
layer10.2 <- geom_text(data = cs.long.merged %>% filter(chromosome == oncogenes[10,]$chr),
                       aes(x=(oncogenes[10,]$middle-18500000), y=0.75), label=oncogenes[10,]$Gene, angle = 90, size = 9)
layer11.2 <- geom_text(data = cs.long.merged %>% filter(chromosome == oncogenes[11,]$chr),
                       aes(x=(oncogenes[11,]$middle+18500000), y=0.60), label=oncogenes[11,]$Gene, angle = 90, size = 9)


layers <- c(layer1, layer2, layer3, layer4, layer5, layer6, layer7, layer8, layer9, layer10, layer11,
            layer1.2, layer2.2, layer3.2, layer4.2, layer5.2, layer6.2, layer7.2, layer8.2, layer9.2, layer10.2, layer11.2)

## Plot them over CNV plot
## For plotting, we select columns 1-7 in order to remove NA values containing columns. Otherwise, cnFreq will remove those rows.
png(
  file = "CNV_genomewide.png",
  res = 200,
  width = 30,
  height = 22,
  units = "in"
)
cnFreq(cs.long.merged[1:7], CN_low_cutoff = 1.5, CN_high_cutoff = 2.5, x_title_size = 20, y_title_size = 20, facet_lab_size = 20, genome = "hg38", plotLayer=layers)
dev.off()

svg(
  file = "CNV_genomewide.svg",
  width = 30,
  height = 22
)
cnFreq(cs.long.merged[1:7], CN_low_cutoff = 1.5, CN_high_cutoff = 2.5, x_title_size = 20, y_title_size = 20, facet_lab_size = 20, genome = "hg38", plotLayer=layers)
dev.off()

## No cell lines 
cs.long.merged.ncl <- filter(cs.long.merged, study != "cell_lines_gabriella_kinker")

png(
  file = "CNV_genomewide_nocelllines.png",
  res = 200,
  width = 30,
  height = 22,
  units = "in"
)
cnFreq(cs.long.merged.ncl[1:7], CN_low_cutoff = 1.5, CN_high_cutoff = 2.5, x_title_size = 20, y_title_size = 20, facet_lab_size = 20, genome = "hg38", plotLayer=layers)
dev.off()


svg(
  file = "CNV_genomewide_nocelllines.svg",
  width = 30,
  height = 22
)
cnFreq(cs.long.merged.ncl[1:7], CN_low_cutoff = 1.5, CN_high_cutoff = 2.5, x_title_size = 20, y_title_size = 20, facet_lab_size = 20, genome = "hg38", plotLayer=layers)

dev.off()

## Only primary tumors
cs.long.merged.p <- filter(cs.long.merged, sample_type == "p", study != "cell_lines_gabriella_kinker")

png(
  file = "CNV_genomewide_onlyprimary.png",
  res = 200,
  width = 45,
  height = 15,
  units = "in"
)
cnFreq(cs.long.merged.p[1:7], CN_low_cutoff = 1.5, CN_high_cutoff = 2.5, x_title_size = 20, y_title_size = 20, facet_lab_size = 20, genome = "hg38", plotLayer=layers)
dev.off()


svg(
  file = "CNV_genomewide_onlyprimary.svg",
  width = 45,
  height = 15
)
cnFreq(cs.long.merged.p[1:7], CN_low_cutoff = 1.5, CN_high_cutoff = 2.5, x_title_size = 20, y_title_size = 20, facet_lab_size = 20, genome = "hg38", plotLayer=layers)

dev.off()

## Only metastatic tumors
cs.long.merged.m <- filter(cs.long.merged, sample_type == "m", study != "cell_lines_gabriella_kinker")

png(
  file = "CNV_genomewide_onlymets.png",
  res = 200,
  width = 45,
  height = 15,
  units = "in"
)
cnFreq(cs.long.merged.m[1:7], CN_low_cutoff = 1.5, CN_high_cutoff = 2.5, x_title_size = 20, y_title_size = 20, facet_lab_size = 20, genome = "hg38", plotLayer=layers)
dev.off()


svg(
  file = "CNV_genomewide_onlymets.svg",
  width = 45,
  height = 15
)
cnFreq(cs.long.merged.m[1:7], CN_low_cutoff = 1.5, CN_high_cutoff = 2.5, x_title_size = 20, y_title_size = 20, facet_lab_size = 20, genome = "hg38", plotLayer=layers)

dev.off()


#############################

## SUPPLEMENTARY FIGURE: By cancer type (only primary, no cell lines) #############
## Get top primary sites by clones
p.table <- as.data.frame(table(cs.long.merged$tumor_type))
p.top <- p.table[order(p.table$Freq, decreasing = T),][1:16,]
top_primary_sites <- as.character(p.top$Var1)

## Test liquid

aml <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "LAML"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
              x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "LAML") + theme(plot.title = element_text(size = 40))

all <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "ALL"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
              x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "ALL") + theme(plot.title = element_text(size = 40))

cll <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "CLL"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
              x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "CLL") + theme(plot.title = element_text(size = 40))

mm <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "MM"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
             x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "MM") + theme(plot.title = element_text(size = 40))

png(
  file = "CNV_genomewide_liquid_tumors.png",
  res = 200,
  width = 45,
  height = 30,
  units = "in"
)

aml + all + cll + mm +
  plot_layout(ncol = 1) & theme(plot.margin = unit(c(0.35,0.05,0.05,0.05), "in"))

dev.off()


#Make CN plots for each cancer type (Incluiding top 10 in fig 1 + other interesting examples)

cnf.laml <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "LAML"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "LAML") + theme(plot.title = element_text(size = 40))

cnf.brca <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "BRCA"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "BRCA") + theme(plot.title = element_text(size = 40))

cnf.luad <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "LUAD"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "LUAD") + theme(plot.title = element_text(size = 40))

cnf.esca <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "ESCA"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "ESCA") + theme(plot.title = element_text(size = 40))

cnf.mm <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "MM"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                 x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "MM") + theme(plot.title = element_text(size = 40))

cnf.gbm <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "GBM"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                  x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "GBM") + theme(plot.title = element_text(size = 40))

cnf.paad <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "PAAD"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD") + theme(plot.title = element_text(size = 40))

cnf.lusc <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "LUSC"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "LUSC") + theme(plot.title = element_text(size = 40))

cnf.coad <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "COAD"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "COAD") + theme(plot.title = element_text(size = 40))

cnf.kirc <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "KIRC"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "KIRC") + theme(plot.title = element_text(size = 40))

cnf.skcm <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "SKCM"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "SKCM") + theme(plot.title = element_text(size = 40))

cnf.nb <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "NB"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                 x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "NB") + theme(plot.title = element_text(size = 40))

png(
  file = "CNV_genomewide_by_tumor_site.png",
  res = 200,
  width = 45,
  height = 60,
  units = "in"
)

cnf.brca + cnf.coad + cnf.esca + cnf.gbm + cnf.kirc + cnf.laml + cnf.luad + cnf.lusc + cnf.mm + cnf.nb + cnf.paad + cnf.skcm +
  plot_layout(ncol = 1) & theme(plot.margin = unit(c(0.35,0.05,0.05,0.05), "in"))

dev.off()

svg(
  file = "CNV_genomewide_by_tumor_site.svg",
  width = 45,
  height = 60
)

cnf.brca + cnf.coad + cnf.esca + cnf.gbm + cnf.kirc + cnf.laml + cnf.luad + cnf.lusc + cnf.mm + cnf.nb + cnf.paad + cnf.skcm +
  plot_layout(ncol = 1) & theme(plot.margin = unit(c(0.35,0.05,0.05,0.05), "in"))

dev.off()


cnf.prad <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "PRAD"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PRAD") + theme(plot.title = element_text(size = 40))

cnf.sarc <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "SARC"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "SARC") + theme(plot.title = element_text(size = 40))

cnf.ov <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "OV"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                 x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "OV") + theme(plot.title = element_text(size = 40))

cnf.ucec <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "UCEC"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "UCEC") + theme(plot.title = element_text(size = 40))

cnf.bcc <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "BCC"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                  x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "BCC") + theme(plot.title = element_text(size = 40))

cnf.hnsc <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "HNSC"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "HNSC") + theme(plot.title = element_text(size = 40))

cnf.sclc <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "SCLC"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "SCLC") + theme(plot.title = element_text(size = 40))

cnf.blca <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "BLCA"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "BLCA") + theme(plot.title = element_text(size = 40))

cnf.read <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "READ"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "READ") + theme(plot.title = element_text(size = 40))

cnf.chol <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "CHOL"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "CHOL") + theme(plot.title = element_text(size = 40))

cnf.escc <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "ESCC"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "ESCC") + theme(plot.title = element_text(size = 40))

cnf.cll <- cnFreq(filter(cs.long.merged.p[1:7], tumor_type == "CLL"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                  x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "CLL") + theme(plot.title = element_text(size = 40))


png(
  file = "CNV_genomewide_by_tumor_site_others.png",
  res = 200,
  width = 45,
  height = 60,
  units = "in"
)

cnf.prad + cnf.sarc + cnf.ov + cnf.ucec + cnf.bcc + cnf.hnsc + cnf.sclc + cnf.blca + cnf.read + cnf.chol + cnf.escc + cnf.cll +
  plot_layout(ncol = 1) & theme(plot.margin = unit(c(0.35,0.05,0.05,0.05), "in"))

dev.off()

##########################################

## FIGURE 2: Table with most recurrent amp and del by TC (no cell lines) ###################
all_m <- list()

for (n in 1:10) {  
  cs.m <- filter(cs.long.merged.ncl, cluster == n)
  all_m[[n]] <- cs.m  
}

top_amp_del <- data.frame(region_name = character(),
                          region = character(),
                          freq.amp = double(),
                          freq.del = double(),
                          type = character())

for(m in all_m){
  number <- unique(m$cluster)
  rnk <- m %>%
    group_by(region_name)
  rnk <- rnk %>% summarise(freq.amp = sum(segmean >= 2.5, na.rm=TRUE)/length(unique(m$sample)),
                           freq.del = sum(segmean <= 1.5, na.rm=TRUE)/length(unique(m$sample)))
  top_amp <- filter(rnk, freq.amp >= 0.5) %>%
    mutate(type = "a")
  top_del <- filter(rnk, freq.del >= 0.5) %>%
    mutate(type = "d")
  
  tops <- rbind(top_amp, top_del) %>%
    mutate(cluster = number)
  
  top_amp_del <- rbind(top_amp_del, tops)
  
}

top_amp_del <- top_amp_del[order(top_amp_del$cluster, top_amp_del$type, top_amp_del$region_name), ]
top_amp_del <- top_amp_del %>%
  mutate(cnv = paste0(type, "_", region_name))
top_amp_del$duplicated <- duplicated(top_amp_del$cnv) | duplicated(top_amp_del$cnv, fromLast = TRUE)
top_amp_del <- mutate(top_amp_del,
                      specific = ifelse(duplicated, "No", "Yes"),
                      arm = sub("(^chr\\d+[pq])\\d+.*", "\\1", region_name))


write.table(top_amp_del, "/local/bc_meta/scevan/cytobands_CNV/top_amp_del_therapeutic_clusters_ncl.tsv", sep = "\t", row.names = FALSE)

ns <- top_amp_del %>%
  filter(specific == "No") %>%
  group_by(cnv) %>%
  mutate(all_clusters = paste(unique(cluster), collapse = ", ")) %>%
  ungroup() %>%
  dplyr::select(cnv, all_clusters) %>%
  unique()


s <- top_amp_del %>%
  filter(specific == "Yes")

c4 <- filter(s, cluster == 4, !(arm %in% c("chr3p", "chr3q", "chr5p", "chr5q", "chr8p", "chr8q"))) #Typical from esophageal
c5 <- filter(s, cluster == 5, !(arm %in% c("chr10p", "chr10q", "chr7p", "chr7q"))) #Typical from glioblastoma

## Make plots 
all_cnf <- list()

for (n in 1:10) {
  cnf <- cnFreq(all_m[[n]][1:7], CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
               x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") +
    labs(title = paste0("Cluster ", n)) + 
    theme(plot.title = element_text(size = 40))
  all_cnf[[n]] <- cnf
}

all_cnf_combined <- wrap_plots(all_cnf) + 
  plot_layout(ncol = 1) & 
  theme(plot.margin = unit(c(0.35, 0.05, 0.05, 0.05), "in"))

png(
  file = "/local/bc_meta/scevan/CNV_genomewide_by_tc.png",
  res = 100,
  width = 30,
  height = 40,
  units = "in"
)

all_cnf_combined

dev.off()


png(
  file = "/local/bc_meta/scevan/CNV_genomewide_by_tc_4and5.png",
  res = 100,
  width = 40,
  height = 15,
  units = "in"
)

wrap_plots(all_cnf[4:5]) + 
  plot_layout(ncol = 1) & 
  theme(plot.margin = unit(c(0.35, 0.05, 0.05, 0.05), "in"))

dev.off()


## General shared CNV (CELL LINES INCLUDED)
rnk <- cs.long.merged %>%
  group_by(region_name) %>% 
  summarise(freq.amp = sum(segmean >= 2.5, na.rm=TRUE)/length(unique(cs.long.merged$sample)),
            freq.del = sum(segmean <= 1.5, na.rm=TRUE)/length(unique(cs.long.merged$sample)))

top_amp <- filter(rnk, freq.amp >= 0.35) %>%
  mutate(type = "a")
top_del <- filter(rnk, freq.del >= 0.35) %>%
  mutate(type = "d")

top_amp_del_general <- rbind(top_amp, top_del)

## Add frequency by sample

freq_df_amp <- cs.long.merged %>%
  filter(region_name %in% filter(top_amp_del_general, type == "a")$region_name) %>%
  filter(segmean >= 2.5) %>%
  group_by(region_name) %>%
  summarise(sample_freq = length(unique(study__sample)) / length(unique(cs.long.merged$study__sample)))

freq_df_del <- cs.long.merged %>%
  filter(region_name %in% filter(top_amp_del_general, type == "d")$region_name) %>%
  filter(segmean <= 1.5) %>%
  group_by(region_name) %>%
  summarise(sample_freq = length(unique(study__sample)) / length(unique(cs.long.merged$study__sample)))

top_amp_del_general_sample <- rbind(freq_df_amp, freq_df_del)


top_amp_del_general <- left_join(top_amp_del_general, top_amp_del_general_sample, by = "region_name")

write.table(top_amp_del_general, "/local/bc_meta/scevan/cytobands_CNV/top_amp_del_global.tsv", sep = "\t", row.names = FALSE)

## Test grouped recurrent amplifications

get_sample_pct <- function(names, type){
  if (type == "a") {
    chr <- cs.long.merged %>%
      group_by(region_name) %>%
      filter(region_name %in% names & segmean >= 2.5)
    
    pct1 <- length(unique(chr$sample))/length(unique(cs.long.merged$sample))
    pct2 <- length(unique(chr$study__sample))/length(unique(cs.long.merged$study__sample))
  }
  
  if (type == "d") {
    chr <- cs.long.merged %>%
      group_by(region_name) %>%
      filter(region_name %in% names & segmean <= 1.5)
    
    pct1 <- length(unique(chr$sample))/length(unique(cs.long.merged$sample))
    pct2 <- length(unique(chr$study__sample))/length(unique(cs.long.merged$study__sample))
  }
  
  cat(paste0("% subclones: ", pct1, "\n% samples: ", pct2))
}


get_sample_pct(names = c("chr10p12.31", "chr10p12.32", "chr10p12.33", "chr10p13"), type = "d")
get_sample_pct(names = c("chr12q13.13", "chr12q13.2"), type = "a")
get_sample_pct(names = c("chr12p13.2", "chr12p13.31"), type = "d")
get_sample_pct(names = c("chr1q21.3", "chr1q22"), type = "a")
get_sample_pct(names = c("chr8p21.3"), type = "d")

## General shared CNV (NO CELL LINES)
rnk <- cs.long.merged.ncl %>%
  group_by(region_name) %>% 
  summarise(freq.amp = sum(segmean >= 2.5, na.rm=TRUE)/length(unique(cs.long.merged.ncl$sample)),
            freq.del = sum(segmean <= 1.5, na.rm=TRUE)/length(unique(cs.long.merged.ncl$sample)))

top_amp <- filter(rnk, freq.amp >= 0.35) %>%
  mutate(type = "a")
top_del <- filter(rnk, freq.del >= 0.35) %>%
  mutate(type = "d")

top_amp_del_general <- rbind(top_amp, top_del)

## Add frequency by sample

freq_df_amp <- cs.long.merged.ncl %>%
  filter(region_name %in% filter(top_amp_del_general, type == "a")$region_name) %>%
  filter(segmean >= 2.5) %>%
  group_by(region_name) %>%
  summarise(sample_freq = length(unique(study__sample)) / length(unique(cs.long.merged.ncl$study__sample)))

freq_df_del <- cs.long.merged.ncl %>%
  filter(region_name %in% filter(top_amp_del_general, type == "d")$region_name) %>%
  filter(segmean <= 1.5) %>%
  group_by(region_name) %>%
  summarise(sample_freq = length(unique(study__sample)) / length(unique(cs.long.merged.ncl$study__sample)))

top_amp_del_general_sample <- rbind(freq_df_amp, freq_df_del)


top_amp_del_general <- left_join(top_amp_del_general, top_amp_del_general_sample, by = "region_name")

write.table(top_amp_del_general, "/local/bc_meta/scevan/cytobands_CNV/top_amp_del_global_ncl.tsv", sep = "\t", row.names = FALSE)

## Test grouped recurrent amplifications

get_sample_pct <- function(names, type){
  if (type == "a") {
    chr <- cs.long.merged.ncl %>%
      group_by(region_name) %>%
      filter(region_name %in% names & segmean >= 2.5)

    pct1 <- length(unique(chr$sample))/length(unique(cs.long.merged.ncl$sample))
    pct2 <- length(unique(chr$study__sample))/length(unique(cs.long.merged.ncl$study__sample))
  }
  
  if (type == "d") {
    chr <- cs.long.merged.ncl %>%
      group_by(region_name) %>%
      filter(region_name %in% names & segmean <= 1.5)
    
    pct1 <- length(unique(chr$sample))/length(unique(cs.long.merged.ncl$sample))
    pct2 <- length(unique(chr$study__sample))/length(unique(cs.long.merged.ncl$study__sample))
  }
  
  cat(paste0("% subclones: ", pct1, "\n% samples: ", pct2))
}


get_sample_pct(names = c("chr10p12.31", "chr10p12.32", "chr10p12.33", "chr10p13"), type = "d")
get_sample_pct(names = c("chr12q13.13", "chr12q13.2"), type = "a")
get_sample_pct(names = c("chr12p13.2", "chr12p13.31"), type = "d")
get_sample_pct(names = c("chr1q21.3", "chr1q22"), type = "a")
get_sample_pct(names = c("chr8p21.3"), type = "d")

##################
## PDAC USECASE

pdac_by_subclone <- cs.long.merged %>%
  filter(tumor_type == "PAAD", sample_type != "") %>% 
  dplyr::select(sample, study, sample_original, cluster, stage, sample_type, treated, refined_tumor_site) %>%
  unique()

pdac_by_sample <- cs.long.merged %>%
  filter(tumor_type == "PAAD", sample_type != "") %>% 
  dplyr::select(study, sample_original, cluster, stage, sample_type, treated, refined_tumor_site) %>%
  unique()

table(pdac_by_subclone$cluster, pdac_by_subclone$sample_type)

table(pdac_by_sample$cluster, pdac_by_sample$sample_type)


cnf.pdac1p <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="1", sample_type =="p"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                   x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 1 (primary)") + theme(plot.title = element_text(size = 40))

cnf.pdac1m <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="1", sample_type =="m"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                 x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 1 (met)") + theme(plot.title = element_text(size = 40))

cnf.pdac2p <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="2", sample_type =="p"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 2 (primary)") + theme(plot.title = element_text(size = 40))

cnf.pdac2m <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="2", sample_type =="m"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 2 (met)") + theme(plot.title = element_text(size = 40))

cnf.pdac3p <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="3", sample_type =="p"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 3 (primary)") + theme(plot.title = element_text(size = 40))

cnf.pdac3m <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="3", sample_type =="m"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 3 (met)") + theme(plot.title = element_text(size = 40))

cnf.pdac4p <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="4", sample_type =="p"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 4 (primary)") + theme(plot.title = element_text(size = 40))

cnf.pdac6p <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="6", sample_type =="p"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 6 (primary)") + theme(plot.title = element_text(size = 40))

cnf.pdac7p <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="7", sample_type =="p"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 7 (primary)") + theme(plot.title = element_text(size = 40))

cnf.pdac7m <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="7", sample_type =="m"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 7 (met)") + theme(plot.title = element_text(size = 40))

cnf.pdac9p <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="9", sample_type =="p"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 9 (primary)") + theme(plot.title = element_text(size = 40))

cnf.pdac9m <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="9", sample_type =="m"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 9 (met)") + theme(plot.title = element_text(size = 40))

cnf.pdac10p <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="10", sample_type =="p"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 10 (primary)") + theme(plot.title = element_text(size = 40))


png(
  file = "CNV_PAAD_usecase2.png",
  res = 200,
  width = 45,
  height = 60,
  units = "in"
)

cnf.pdac1p + cnf.pdac1m + cnf.pdac2p + cnf.pdac2m + cnf.pdac3p + cnf.pdac3m + cnf.pdac4p + cnf.pdac6p + cnf.pdac7p + cnf.pdac7m + 
  cnf.pdac9p + cnf.pdac9m + cnf.pdac10p +
  plot_layout(ncol = 1) & theme(plot.margin = unit(c(0.35,0.05,0.05,0.05), "in"))

dev.off()


cnf.pdac1 <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="1"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                     x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 1") + theme(plot.title = element_text(size = 40))

cnf.pdac2 <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="2"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                    x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 2") + theme(plot.title = element_text(size = 40))

cnf.pdac3 <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="3"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                    x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 3") + theme(plot.title = element_text(size = 40))

cnf.pdac4 <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="4"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                    x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 4") + theme(plot.title = element_text(size = 40))

cnf.pdac6 <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="6"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                    x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 6") + theme(plot.title = element_text(size = 40))

cnf.pdac7 <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="7"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                    x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 7") + theme(plot.title = element_text(size = 40))

cnf.pdac9 <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="9"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                    x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 9") + theme(plot.title = element_text(size = 40))

cnf.pdac10 <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster =="10"), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                    x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 10") + theme(plot.title = element_text(size = 40))




png(
  file = "CNV_PAAD_usecase3.png",
  res = 200,
  width = 45,
  height = 40,
  units = "in"
)

cnf.pdac1 + cnf.pdac2 + cnf.pdac3 + cnf.pdac4 + cnf.pdac6 + cnf.pdac7 + cnf.pdac9 + cnf.pdac10 +
  plot_layout(ncol = 1) & theme(plot.margin = unit(c(0.35,0.05,0.05,0.05), "in"))

dev.off()


cnf.pdac123 <- cnFreq(filter(cs.long.merged, tumor_type == "PAAD", cluster %in% c("1", "2", "3")), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                    x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TCs 1, 2 and 3") + theme(plot.title = element_text(size = 40))


png(
  file = "CNV_PAAD_usecase4.png",
  res = 200,
  width = 45,
  height = 20,
  units = "in"
)

cnf.pdac123 + cnf.pdac7 +
  plot_layout(ncol = 1) & theme(plot.margin = unit(c(0.35,0.05,0.05,0.05), "in"))

dev.off()



paad.123.p <- filter(cs.long.merged, tumor_type == "PAAD", cluster %in% c("1", "2", "3"), sample_type == "p")
paad.123.p$group <- "123p"
paad.123.m <- filter(cs.long.merged, tumor_type == "PAAD", cluster %in% c("1", "2", "3"), sample_type == "m")
paad.123.m$group <- "123m"
paad.7.p <- filter(cs.long.merged, tumor_type == "PAAD", cluster =="7", sample_type =="p")
paad.7.p$group <- "7p"
paad.7.m <- filter(cs.long.merged, tumor_type == "PAAD", cluster =="7", sample_type =="m")
paad.7.m$group <- "7m"


cnf.pdac123p <- cnFreq(paad.123.p, CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                      x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TCs 1, 2 and 3 (primary)") + theme(plot.title = element_text(size = 40))

cnf.pdac123m <- cnFreq(paad.123.m, CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                       x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TCs 1, 2 and 3 (met)") + theme(plot.title = element_text(size = 40))

cnf.pdac7p <- cnFreq(paad.7.p, CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                       x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TC 7 (primary)") + theme(plot.title = element_text(size = 40))

cnf.pdac7m <- cnFreq(paad.7.m, CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
                       x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "PAAD in TCs 7 (met)") + theme(plot.title = element_text(size = 40))


png(
  file = "CNV_PAAD_usecase4.png",
  res = 200,
  width = 30,
  height = 20,
  units = "in"
)

cnf.pdac123p + cnf.pdac123m + cnf.pdac7p + cnf.pdac7m +
  plot_layout(ncol = 1) & theme(plot.margin = unit(c(0.35,0.05,0.05,0.05), "in"))

dev.off()



all_m <- list(paad.123.p, paad.123.m, paad.7.p, paad.7.m)

top_amp_del <- data.frame(region_name = character(),
                          region = character(),
                          freq.amp = double(),
                          freq.del = double(),
                          type = character())


for(m in all_m){
  number <- unique(m$group)
  rnk <- m %>%
    group_by(region_name)
  rnk <- rnk %>% summarise(freq.amp = sum(segmean >= 2.5, na.rm=TRUE)/length(unique(m$sample)),
                           freq.del = sum(segmean <= 1.5, na.rm=TRUE)/length(unique(m$sample)))
  top_amp <- filter(rnk, freq.amp >= 0.5) %>%
    mutate(type = "a")
  top_del <- filter(rnk, freq.del >= 0.5) %>%
    mutate(type = "d")
  
  tops <- rbind(top_amp, top_del) %>%
    mutate(group = number)
  
  top_amp_del <- rbind(top_amp_del, tops)
  
}

top_amp_del <- top_amp_del[order(top_amp_del$group, top_amp_del$type, top_amp_del$region_name), ]
top_amp_del <- top_amp_del %>%
  mutate(cnv = paste0(type, "_", region_name))
top_amp_del$duplicated <- duplicated(top_amp_del$cnv) | duplicated(top_amp_del$cnv, fromLast = TRUE)
top_amp_del <- mutate(top_amp_del,
                      specific = ifelse(duplicated, "No", "Yes"),
                      arm = sub("(^chr\\d+[pq])\\d+.*", "\\1", region_name))


write.table(top_amp_del, "/local/bc_meta/scevan/cytobands_CNV/top_amp_del_PDAC_usecase.tsv", sep = "\t", row.names = FALSE)

cnv123p <- filter(top_amp_del, group == "123p") %>% dplyr::select(cnv) %>% unlist()
cnv123m <- filter(top_amp_del, group == "123m") %>% dplyr::select(cnv) %>% unlist()
cnv123common <- intersect(cnv123p, cnv123m)


cnv7p <- filter(top_amp_del, group == "7p") %>% dplyr::select(cnv) %>% unlist()
cnv7m <- filter(top_amp_del, group == "7m") %>% dplyr::select(cnv) %>% unlist()
cnv7common <- intersect(cnv7p, cnv7m)

cnv_specific_123m <- setdiff(cnv123m, cnv7m)
cnv_specific_123m <- setdiff(cnv_specific_123m, cnv123p)

cnv_specific_7m <- setdiff(cnv7m, cnv123m)
cnv_specific_7m <- setdiff(cnv_specific_7m, cnv7p)

## Create contingency matrix
contingency_table <- table(top_amp_del$group, top_amp_del$cnv)

groups <- unique(top_amp_del$group)
result_matrix <- matrix(0, nrow = length(groups), ncol = length(groups),
                        dimnames = list(groups, groups))

for (i in 1:length(groups)) {
  for (j in 1:length(groups)) {
    common_cnv <- sum(contingency_table[i, ] > 0 & contingency_table[j, ] > 0)
    result_matrix[i, j] <- common_cnv
  }
}

print(result_matrix)


##################
## BRCA USECASE
cs.long.merged.tc <- filter(cs.long.merged, !(is.na(cluster)))

## BRCA TC10 vs rest
tmp <- cnFreq(filter(cs.long.merged.tc, tumor_type == "BRCA", cluster == 10), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
              x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "TC10 BRCA") + theme(plot.title = element_text(size = 40))

tmp2 <- cnFreq(filter(cs.long.merged.tc, tumor_type == "BRCA", cluster != 10), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5,
               x_title_size = 0, y_title_size = 0, facet_lab_size = 20, genome = "hg38") + labs(title = "Other TCs BRCA") + theme(plot.title = element_text(size = 40))


svg(
  file = "CNV_BRCA_TC10.svg",
  width = 45,
  height = 30
)
tmp + tmp2 +
  plot_layout(ncol = 1) & theme(plot.margin = unit(c(0.35,0.05,0.05,0.05), "in"))

dev.off()


brca10 <- filter(cs.long.merged.tc, tumor_type == "BRCA", cluster == 10)
brca10$group <- "TC10"
brcaTC <- filter(cs.long.merged.tc, tumor_type == "BRCA", cluster != 10)
brcaTC$group <- "Other"


all_m <- list(brca10, brcaTC)

top_amp_del <- data.frame(region_name = character(),
                          region = character(),
                          freq.amp = double(),
                          freq.del = double(),
                          type = character())

for(m in all_m){
  number <- unique(m$group)
  rnk <- m %>%
    group_by(region_name)
  rnk <- rnk %>% summarise(freq.amp = sum(segmean >= 2.5, na.rm=TRUE)/length(unique(m$sample)),
                           freq.del = sum(segmean <= 1.5, na.rm=TRUE)/length(unique(m$sample)))
  top_amp <- filter(rnk, freq.amp > 0.5) %>% # I changed this filter to select tops only (removed the =, as many cnv were on the edge)
    mutate(type = "a")
  top_del <- filter(rnk, freq.del > 0.5) %>%
    mutate(type = "d")
  
  tops <- rbind(top_amp, top_del) %>%
    mutate(group = number)
  
  top_amp_del <- rbind(top_amp_del, tops)
  
}

top_amp_del <- top_amp_del[order(top_amp_del$group, top_amp_del$type, top_amp_del$region_name), ]
top_amp_del <- top_amp_del %>%
  mutate(cnv = paste0(type, "_", region_name))
top_amp_del$duplicated <- duplicated(top_amp_del$cnv) | duplicated(top_amp_del$cnv, fromLast = TRUE)
top_amp_del <- mutate(top_amp_del,
                      specific = ifelse(duplicated, "No", "Yes"),
                      arm = sub("(^chr\\d+[pq])\\d+.*", "\\1", region_name))


write.table(top_amp_del, "/local/bc_meta/scevan/cytobands_CNV/top_amp_del_BRCA_usecase.tsv", sep = "\t", row.names = FALSE)

cnv10 <- filter(top_amp_del, group == "TC10") %>% dplyr::select(cnv) %>% unlist()
cnvTC <- filter(top_amp_del, group == "Other") %>% dplyr::select(cnv) %>% unlist()
common_cnv <- intersect(cnv10, cnvTC)
specific10 <- setdiff(cnv10, cnvTC)
specificTC <- setdiff(cnvTC, cnv10)

specific10_df <-  filter(top_amp_del, cnv %in% specific10)
write.table(specific10_df, "/local/bc_meta/scevan/cytobands_CNV/top_amp_del_BRCA_usecase_tc10_specific.tsv", sep = "\t", row.names = FALSE)

## These are actually:
## AMPS:1p33-34.2, 1q23.3, 8q22.1-24.3, chr8q24.11
## DELS: 13q13-14, 14q21-24, 5q12.3-15, 5q21.1-22.2

## TC4 USECASE

cnv_start <- geom_vline(data = cs.long.merged %>% filter(cluster == 4, chromosome == "chr3"),
                      aes(xintercept = 149200001), col = "#3a3b3c", linetype = 2, size = 1)
cnv_end <- geom_vline(data = cs.long.merged %>% filter(cluster == 4, chromosome == "chr3"),
                      aes(xintercept = 198295559), col = "#3a3b3c", linetype = 2, size = 1)
cnv_name <- geom_text(data = cs.long.merged %>% filter(cluster == 4, chromosome == "chr3"),
                       aes(x=((149200001+198295559)/2), y=0.90), label="3q25-q29", angle = 90, size = 9)

layers <- c(cnv_start, cnv_end, cnv_name)

## Plot them over CNV plot
png(
  file = "CNV_TC4.png",
  res = 200,
  width = 42,
  height = 22,
  units = "in"
)
cnFreq(cs.long.merged %>% filter(cluster == 4), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5, x_title_size = 20, y_title_size = 20, facet_lab_size = 20, genome = "hg38", plotLayer=layers)

dev.off()



svg(
  file = "CNV_TC4.svg",
  width = 35,
  height = 22
)

cnFreq(cs.long.merged %>% filter(cluster == 4), CN_low_cutoff = 1.5, CN_high_cutoff = 2.5, x_title_size = 20, y_title_size = 20, facet_lab_size = 20, genome = "hg38", plotLayer=layers)

dev.off()
