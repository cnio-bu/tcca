library(GenomicRanges)
library(BPCells)
library(tidyverse)
library(plyranges)

setwd("/local/bc_meta/scevan/")

## Load data

mat <- open_matrix_dir("/local/bc_meta/scevan/cnv_segments_clones_lvl2_fullbpcellsmatrix/")

## Get dataframe to build GRanges object
regions <- rownames(mat)
split_reg <- strsplit(regions, "_")
regions_df <- as.data.frame(do.call(rbind, split_reg))
colnames(regions_df) <- c("Chr", "Start", "End")
regions_df <- mutate(regions_df,
                     Chr = paste0("chr", Chr))

gr_scevan <- GRanges(
  seqnames = regions_df$Chr,
  ranges = IRanges(start = as.numeric(regions_df$Start), end = as.numeric(regions_df$End), names = rownames(regions_df)),
  scevan_name =  regions
)

mean(width(gr_scevan)) # Average width is ~20.31 Mb, way bigger than cytobands (~2.5 Mb)


## Make a custom ranges set of the human genome
# Load data of chromosome lengths of hg38 from Genome Reference Consortium
chr_lengths <- read.table("/local/bc_meta/SCEVAN_information/seurat5/chr_lengths.tsv", header = T, sep = "\t")

# Generate the ranges
avg_length <- 20306177

gr_custom <- GRanges()

for (chr in 1:24){
  # Define the total range
  total_range <- 1:chr_lengths[chr,]$bp
  
  # Calculate the number of ranges
  num_ranges <- ceiling(length(total_range)/avg_length)
  
  # Calculate the actual length of each range
  actual_length <- ceiling(length(total_range)/num_ranges)
  
  # Generate the ranges
  ranges <- sapply(1:num_ranges, function(i) {
    start_index <- (i - 1) * actual_length + 1
    end_index <- min(i * actual_length, length(total_range))
    return(c(start_index, end_index))
  })
  
  rownames(ranges) <- c("start", "end")
  ranges <- as.data.frame(t(ranges))
  
  gr_ranges <- GRanges(
    seqnames = Rle(chr_lengths[chr,]$Chromosome, num_ranges),
    ranges = IRanges(start = as.numeric(ranges$start), end = as.numeric(ranges$end)),
    custom_name =  paste0(chr_lengths[chr,]$Chromosome, head(letters, num_ranges))
  )
  
  gr_custom <- c(gr_custom, gr_ranges)
}

## Use plyranges to join both 

joined_gr <- join_overlap_intersect(gr_scevan, gr_custom)
widths <- as.data.frame(ranges(joined_gr))
joined_gr$width <- widths$width

result_gr <- joined_gr %>%
  group_by(scevan_name) %>%
  filter(width == max(width, na.rm = TRUE))


renaming <- data.frame(scevan_name = result_gr$scevan_name,
                       custom_name = result_gr$custom_name)

write.table(renaming, "scevan_to_custom_segments.tsv", sep = "\t")


##Same stuff but with cytobands (much smaller)

library(AnnotationHub)

proxy <- "lserranor@cnio.es"
proxy <- httr::use_proxy(Sys.getenv('http_proxy'))
httr::set_config(proxy)
AnnotationHub::setAnnotationHubOption("PROXY", proxy)
AnnotationHub::getAnnotationHubOption("LOCAL")

hub <- AnnotationHub()

hub_hg38 <- subset(hub, 
                   (hub$species == "Homo sapiens") & (hub$genome == "hg38"))

cytobands  <- hub_hg38[[797]]

mean(width(cytobands)) #~2.5 Mb


## Use plyranges to join both 

joined_gr <- join_overlap_intersect(gr_scevan, cytobands)
widths <- as.data.frame(ranges(joined_gr))
joined_gr$width <- widths$width

result_gr <- joined_gr %>%
  group_by(scevan_name) 
## filter(width == max(width, na.rm = TRUE))This is removed in the case of cytobands, otherwise we have many empty spaces


renaming <- data.frame(scevan_name = result_gr$scevan_name,
                       custom_name = paste0(seqnames(result_gr), result_gr$name))

write.table(renaming, "scevan_to_cytobands.tsv", sep = "\t")