library(openxlsx)
library(biomaRt)
library(dplyr)

setwd("/home/lmgonzalezb/Documents/bc-meta/SCellBow/GBM/beyondcell_target/TMZ_signature/")

# Paper 1:  Ntafoulis I, Kleijn A, Ju J, et al. Ex vivo drug sensitivity screening 
# predicts response to temozolomide in glioblastoma patients and identifies 
# candidate biomarkers. Br J Cancer. 2023;129(8):1327-1338. doi:10.1038/s41416-023-02402-y
# (Figure 5)

# Up-regulated genes
genes_up <- c(
  "PLTP", "NBPF11", "LRCH3", "FTX",
  "CCDC191", "PWAR5", "PTOV1-AS1", "MAU2",
  "ZNF251", "CCDC30", "ZNF81", "CRTC1",
  "HSPA6", "PPP1R42", "SUDS3", "LENG8",
  "CLN5", "CLEC2D", "PIP4K2B", "ARHGAP31",
  "THAP2", "ARHGAP39", "TSNAXIP1", "SRGAP2",
  "SPPL2B", "MVB12B", "TAOK3", "MTERF2",
  "VMAC", "PCSK5", "JUND", "PTX3",
  "RAP2B", "PCMTD2", "CCDC130", "PWAR6",
  "MIR600HG", "DISC1", "KCTD16", "KIF27",
  "PTCH2", "SOCS7", "PIGM", "SLC16A2",
  "AMY2B", "C4orf47", "KANSL1L", "ADCY10P1",
  "BBS1", "CFAP44"
)

# Down-regulated genes
genes_down <- c(
  "CCNB2", "PFKP", "DMC1", "RPL4P5",
  "MGMT", "RFWD3", "MTFMT", "KIF23",
  "CCZ1", "COX5A", "CENPN", "ZWILCH",
  "CDR2", "ANLN", "SEC11A", "TIPIN",
  "EIF1AXP1", "RPS17", "CIAPIN1", "SNRNP27",
  "PDAP1", "RACGAP1P", "OR7E38P", "CHEK1",
  "NCAPG2", "COMMD4", "PRC1", "RPL4",
  "NSUN5", "MZT2A", "ANP32B", "SNHG16",
  "HAUS4", "GSE1", "HMGA1", "FN3KRP",
  "EIF4H", "RCC1L", "ARNTL2", "AGK",
  "PM20D2", "CHCHD3", "CKAP2L", "SLC35F2",
  "CCNE2", "TUB8", "SETD3", "RRM2",
  "MRM2", "TPMTP1"
)

signatures <- c()
sig1_up <- paste(c("TEMOZOLOMIDE_SENSITIVE_NTAFOULIS_UP",
                   "Ioannis_Ntafoulis_2023",
                   rev(up_genes)), 
                 collapse = "\t")
sig1_down <- paste(c("TEMOZOLOMIDE_SENSITIVE_NTAFOULIS_DOWN",
                     "Ioannis_Ntafoulis_2023",
                     rev(down_genes)), 
                   collapse = "\t")
signatures <- c(signatures, sig1_up, sig1_down)

## Paper 2: Cai HQ, Liu AS, Zhang MJ, et al. Identifying Predictive Gene Expression 
# and Signature Related to Temozolomide Sensitivity of Glioblastomas. Front Oncol. 
# 2020;10:669. Published 2020 May 22. doi:10.3389/fonc.2020.00669 (Suplemmentary
# table 1 and 2).
tab1 <- read.xlsx("Table 1.xlsx", sheet = "sheet1")
# Remove version numbers (keep only the part before the dot)
tab1$ensg_ids <- sub("\\..*", "", tab1$Gene)

# biomaRt to convert ENSG to Gene Symbol
listEnsembl()
ensembl <- useEnsembl((biomart = "genes"))
datasets <- listDatasets(ensembl)

ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attr <- listAttributes((ensembl.con))
filters <- listFilters(ensembl.con)

ensg_gensymbol <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                        filters = "ensembl_gene_id",
                        values = tab1$ensg_ids,
                        mart = ensembl.con)
  
tab1 <- tab1 %>%
  left_join(ensg_gensymbol, by = c("ensg_ids" = "ensembl_gene_id"))

genes_up <- tab1 %>%
  filter(Coefficient >= 0.4) %>%
  filter(!(is.na(external_gene_name) | external_gene_name == "")) %>%
  pull(external_gene_name)

genes_down <- tab1 %>%
  filter(Coefficient <= -0.4) %>%
  filter(!(is.na(external_gene_name) | external_gene_name == "")) %>%
  pull(external_gene_name)

# Add signature 2
sig2_up <- paste(c("TEMOZOLOMIDE_SENSITIVE_CAIHQ_UP",
                   "Hong-Qing_Cai_2020",
                   genes_up), 
                 collapse = "\t")
sig2_down <- paste(c("TEMOZOLOMIDE_SENSITIVE_CAIHQ_DOWN",
                     "Hong-Qing_Cai_2020",
                     genes_down), 
                   collapse = "\t")
signatures <- c(signatures, sig2_up, sig2_down)


# Create a signature only with the protein-coding genes from Cai HQ, et al, 2020
tab2 <- read.xlsx("Table 2.xlsx", sheet = "sheet1")
tab2 <- tab2 %>%
  select(-patient.number.of.high.expression) %>%
  left_join(select(tab1, Gene, external_gene_name), by = c("Gene" = "Gene"))

genes_up <- tab2 %>%
  filter(HR > 1) %>%
  filter(!(is.na(external_gene_name) | external_gene_name == "")) %>%
  pull(external_gene_name)

genes_down <- tab2 %>%
  filter(HR < 1) %>%
  filter(!(is.na(external_gene_name) | external_gene_name == "")) %>%
  pull(external_gene_name)

# Add signature 2
sig3_up <- paste(c("TEMOZOLOMIDE_SENSITIVE_CAIHQ_CODING_GENES_UP",
                   "Hong-Qing_Cai_2020",
                   genes_up), 
                 collapse = "\t")
sig3_down <- paste(c("TEMOZOLOMIDE_SENSITIVE_CAIHQ_CODING_GENES_DOWN",
                     "Hong-Qing_Cai_2020",
                     genes_down), 
                   collapse = "\t")
signatures <- c(signatures, sig3_up, sig3_down)
signatures_lines <- paste(signatures, collapse = "\n")

writeLines(text = signatures_lines, con = "TMZ_signatures.gmt")

  
  

