library(dplyr)
library(tidyverse)

setwd("/storage/scratch01/shared/projects/bc-meta/")

## Generate metadata + clonality info for the full cohort
tcca_metadata <- read.table(
  "single_cell/seurat/tcca/tcca_annotation_raw.tsv",
  header = TRUE,
  sep = "\t"
)
colnames(tcca_metadata)[1] <- "original_barcode"
tcca_metadata <- tcca_metadata %>%
  mutate(malignancy = ifelse(malignancy == "True", TRUE, FALSE))

clonality <- read.table("./single_cell/cna_metadata/full_clonality_table_lvl2.tsv",
                        row.names = NULL)

duplicated <- clonality %>%
  filter(original_barcode %in% original_barcode[duplicated(original_barcode)])

clonality <- clonality %>%
  filter(!(original_barcode %in% duplicated$original_barcode)) %>% # Remove 16 cells with duplicated names
  dplyr::select(original_barcode,
                scevan_prediction,
                confidentNormal,
                scevan_subclone)

full_metadata <- left_join(tcca_metadata, clonality, by = "original_barcode")

## Keep only cells in agreement for better quality of the prediction
full_metadata <- full_metadata %>%
  filter(
    malignancy == FALSE &
      scevan_prediction == "normal" |
      malignancy == TRUE & scevan_prediction == "tumor"
  )

full_metadata <- full_metadata %>%
  mutate(
    subclone_name = paste(study, sample, scevan_subclone, sep = "."),
    ensemble_output = ifelse(malignancy == FALSE, "healthy", subclone_name)
  ) %>%
  column_to_rownames(var = "original_barcode")

## Keep only subclones (malignant cells)
full_metadata <- full_metadata %>%
  filter(malignancy == TRUE)

# Keep metadata at subclonal level
full_metadata <- full_metadata %>%
  select(sample,
         colnames(full_metadata)[8:27],
         tme_archetype,
         subclone_name) %>%
  distinct()

# Add broad cancer type annotation
cancer_type <- list(
  "Brain Cancer" = c("GBM", "MB", "OGD"),
  "Neuroblastic Tumors" = c("GNB", "NB"),
  "Blood Cancer" = c("ALL", "LAML", "CLL", "MM"),
  "Skin Cancer" = c("BCC", "SKCM", "SKSC", "SKAM", "UVM"),
  "Sarcoma/Soft Tissue Cancer" = c("SARC", "GIST", "MESO"),
  "Breast Cancer" = c("BRCA"),
  "Lung Cancer" = c("SCLC", "NSCLC", "LUAD", "LUSC", "LCLC", "PLEU"),
  "Ovarian Cancer" = c("OV"),
  "Colon/Colorectal Cancer" = c("COAD", "READ"),
  "Endometrial/Uterine Cancer" = c("CESC", "UCEC", "UCS"),
  "Liver/Biliary Cancer" = c("LIHC", "CHOL"),
  "Bladder Cancer" = c("BLCA"),
  "Head and Neck Cancer" = c("HNSC"),
  "Prostate Cancer" = c("PRAD"),
  "Kidney Cancer" = c("KRCC", "KTCC", "KIRC", "KIRCH"),
  "Esophageal Cancer" = c("ESCA", "ESCC"),
  "Pancreatic Cancer" = c("PAAD"),
  "Thyroid Cancer" = c("THCA"),
  "Gastric Cancer" = c("STAD"),
  "Miscellaneous Cancer" = c("MISC")
)

cancer_type <- enframe(cancer_type, name = "broad_tumor_type", value = "tumor_type") %>%
  unnest()

full_metadata <-  full_metadata %>%
  left_join(cancer_type, by = "tumor_type") %>%
  relocate(broad_tumor_type, .after = age)

# Add drugs
drug_prediction <- read.table("./single_cell/sctherapy/full_table_drug_prediction.tsv")
drug_prediction <- drug_prediction %>%
  select(-c(cid, Sample))

subclone_metadata <- full_metadata %>%
  right_join(drug_prediction, by = c("subclone_name" = "Subclone")) %>%
  relocate(subclone_name, .before = 1)

# Match MoAs with our custom MoAs
MoAs <- read.table("./reference/final_moas - Collapsed.tsv",
                   header = TRUE,
                   sep = "\t") %>%
  select(preferred.drug.names, collapsed.MoAs) %>%
  distinct()

MoAs_sctherapy <- subclone_metadata %>%
  select(Drug_Name, MoA) %>%
  distinct() %>%
  mutate(Drug_Name = toupper(Drug_Name)) %>%
  column_to_rownames(var = "Drug_Name")

subclone_metadata <- subclone_metadata %>%
  mutate(Drug_Name = toupper(Drug_Name)) %>%
  left_join(MoAs, by = c("Drug_Name" = "preferred.drug.names"))

drug_na <- unique(subclone_metadata$Drug_Name[is.na(subclone_metadata$collapsed.MoAs)])

subclone_metadata <- subclone_metadata %>%
  mutate(collapsed.MoAs = ifelse(is.na(collapsed.MoAs), MoAs_sctherapy[Drug_Name, "MoA"], collapsed.MoAs))

# Adapt MoA for drugs with scTherapy MoA
subclone_metadata <- subclone_metadata %>%
  mutate(
    collapsed.MoAs = case_when(
      collapsed.MoAs == "thymidylate_synthase_inhibitor" ~ "DNA related agent",
      collapsed.MoAs == "tubulin_polymerization_inhibitor" ~ "Microtubule agent",
      collapsed.MoAs == "MDM_inhibitor" ~ "p53 activator/MDM2 inhibitor",
      collapsed.MoAs == "ALK_tyrosine_kinase_receptor_inhibitor" ~ "Kinase inhibitor",
      collapsed.MoAs == "Brutons_tyrosine_kinase_BTK_inhibitor" ~ "Kinase inhibitor",
      collapsed.MoAs == "inosine_monophosphate_dehydrogenase_inhibitor" ~ "DNA related agent",
      collapsed.MoAs == "RAF_inhibitor" ~ "BRAF inhibitor",
      collapsed.MoAs == "RAF_inhibitor" ~ "Cell cycle arrest",
      collapsed.MoAs == "BRD-K24576554" ~ "Other",
      collapsed.MoAs == "BRD-K60237333" ~ "Other",
      collapsed.MoAs == "BRD-K95142244" ~ "Other",
      collapsed.MoAs == "MAP_kinase_inhibitor" ~ "MAPK inhibitor",
      collapsed.MoAs == "adenosine_receptor_antagonist" ~ "Other",
      collapsed.MoAs == "centromere_associated_protein_inhibitor" ~ "Microtubule agent",
      collapsed.MoAs == "PLK_inhibitor" ~ "Cell cycle arrest",
      collapsed.MoAs == "HDAC_inhibitor" ~ "Chromatin agent",
      collapsed.MoAs == "CDK_inhibitor" ~ "Cell cycle arrest",
      collapsed.MoAs == "topoisomerase_inhibitor" ~ "DNA related agent",
      collapsed.MoAs == "DNA_synthesis_inhibitor_microtubule_inhibitor" ~ "Microtubule agent",
      collapsed.MoAs == "kinesin_inhibitor_kinesinlike_spindle_protein_inhibitor" ~ "Microtubule agent",
      collapsed.MoAs == "VEGFR_inhibitor" ~ "VEGFR inhibitor",
      collapsed.MoAs == "potassium_channel_activator" ~ "Other",
      collapsed.MoAs == "KIT_inhibitor_PDGFR_tyrosine_kinase_receptor_inhibitor_VEGFR_inhibitor" ~ "VEGFR inhibitor",
      collapsed.MoAs == "dehydrogenase_inhibitor_inositol_monophosphatase_inhibitor" ~ "Other",
      collapsed.MoAs == "FLT3_inhibitor_JAK_inhibitor" ~ "JAK-STAT signaling inhibitor",
      collapsed.MoAs == "Ttype_calcium_channel_blocker" ~ "Other",
      collapsed.MoAs == "BcrAbl_kinase_inhibitor_FLT3_inhibitor_PDGFR_tyrosine_kinase_receptor_inhibitor" ~ "BCR-ABL inhibitor",
      collapsed.MoAs == "estrogen_receptor_antagonist_selective_estrogen_receptor_modulator_SERM" ~ "Other",
      collapsed.MoAs == "retinoid_receptor_agonist" ~ "Other",
      collapsed.MoAs == "HSP_inhibitor" ~ "HSP inhibitor",
      collapsed.MoAs == "RNA_polymerase_inhibitor" ~ "Transcription inhibitor",
      collapsed.MoAs == "PI3K_inhibitor" ~ "PI3K/AKT/mTOR signaling inhibitor",
      collapsed.MoAs == "JAK_inhibitor" ~ "JAK-STAT signaling inhibitor",
      collapsed.MoAs == "mTOR_inhibitor_PI3K_inhibitor" ~ "PI3K/AKT/mTOR signaling inhibitor",
      collapsed.MoAs == "Aurora_kinase_inhibitor_FLT3_inhibitor_VEGFR_inhibitor" ~ "Other",
      collapsed.MoAs == "BRAF inhibitor;VEGFR inhibitor" ~ "Kinase inhibitor",
      collapsed.MoAs == "BCR-ABL inhibitor;SRC inhibitor" ~ "BCR-ABL inhibitor",
      collapsed.MoAs == "VEGFR inhibitor;MET inhibitor" ~ "Kinase inhibitor",
      collapsed.MoAs == "CDK_inhibitor_cell_cycle_inhibitor_MCL1_inhibitor" ~ "Cell cycle arrest",
      Drug_Name == "CYT387" ~ "JAK-STAT signaling inhibitor",
      collapsed.MoAs == "-" ~ "Other",
      TRUE ~ collapsed.MoAs
    )
  ) %>%
  select(-MoA)


# Add clusters
sctherapy_cluster <- readRDS("./single_cell/sctherapy/results/speclustering_reordered.rds")
sctherapy_cluster_df <- as.data.frame(sctherapy_cluster) %>%
  rownames_to_column(var = "subclone_name")

subclone_metadata <- subclone_metadata %>%
  left_join(sctherapy_cluster_df, by = "subclone_name") %>%
  relocate(Drug_Name,
           Response,
           collapsed.MoAs,
           Dose.nM.,
           Toxicity,
           sctherapy_cluster,
           .after = 1)

colnames(subclone_metadata) <- c(
  "Subclone Name",
  "Drug Name",
  "Response",
  "Drug Mechanism of Action",
  "Dose.nM.",
  "Toxicity",
  "scTherapy Cluster",
  "Sample",
  "Study",
  "Patient",
  "Sex",
  "Age",
  "Broad Tumor Type",
  "Tumor Type",
  "Tumor Subtype",
  "Stage",
  "TNM",
  "Tumor Site",
  "Sample Type",
  "Treated",
  "Treatment Type",
  "Treatment Response",
  "Treatment Info",
  "OS",
  "PFI",
  "Sequencing Tech",
  "Genome Assembly",
  "Data PMID",
  "Refined Tumor Site",
  "TME Archetype"
)

write.table(
  subclone_metadata,
  "./single_cell/sctherapy/results/subclone_level_annotated.tsv",
  sep = "\t",
  row.names = FALSE
)

subclone_final <- subclone %>%
  left_join(select(mps, c("scevan_subclone", "top_MP_clean")), by = c("Subclone.Name" = "scevan_subclone"))

subclone_final <- subclone_final %>% 
  mutate(
    Functional_metaprogram_family = str_extract(Functional_metaprogram, "(?<=_)[^.]+"),
    Functional_metaprogram_family = case_when(
    str_detect(Functional_metaprogram, "LineageSpecific.Hemato") ~ "LineageSpecific.Hemato",
    str_detect(Functional_metaprogram, "LineageSpecific.Neural") ~ "LineageSpecific.Neural",
    str_detect(Functional_metaprogram, "Secretory|Melanocyte") ~ "LineageSpecific.Other",
    TRUE ~ Functional_metaprogram_family
  )
  )
