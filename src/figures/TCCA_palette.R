## TCCA COLOR PALETTES

# Tumor type x15 (14 + Other)
tumor_sites_colors <- c("Bone marrow" = "#FB6467",
                        "Brain" = "#D3C3E0",
                        "Adrenal gland" = "#FF9A4A",
                        "Breast" = "#FF8FAB",
                        "Skin" = "#CA9A8C",
                        "Lung" = "#A2F0EF",
                        "Soft tissue" = "#FFE072",
                        "Esophagus" = "#77C9ED", 
                        "Bladder" = "#58B368",
                        "Lymph node" = "#B6F884",
                        "Liver" = "#309975",
                        "Pancreas" = "#B47EB3",
                        "Ovary" = "#EC8EED",
                        "Prostate" = "#5e91ab",
                        "Colon" = "#BD081A",
                        "Kidney"  = "#DBC586",
                        "Other" = "#6b6363")

# Patient/cell line
patient_ccl_colors <- c("Patient" = "#427394", 
                        "Cell line" = "#F78C1F")

# Treated/Untreated 
treatment_colors <- c("Treated" = "#6ED1BC",
                      "Untreated" = "#D18B6E")

# Adult/Pediatric 
age_colors <- c("Adult" = "#7689DE",
                "Pediatric" = "#a9dce3")

# Male/Female 
sex_colors <- c("Male" = "#ec5c44",
               "Female" = "#1c8c8c")

# Solid/Liquid 
sl_colors <- c("Solid" = "#F6Bd60",
                "Liquid" = "#706695")

# Primary/Met 
pm_colors <- c("Metastasis" = "#C10044",
               "Primary" = "#F0BFD0")

# Therapeutic clusters x5
tcs_colors <- c("5" = "#A3A500", 
                "4"= "#F8766D",
                "3" = "#FFD64C",
                "2" = "#00BF7D",
                "1" = "#00B0F6",
                "0" = "#E76BF3")

# Expression clusters x10
ecs_colors <- c("10" ="#A7D676",
                "9" = "#43978D",
                "8" = "#90E2ED",
                "7" = "#369CBB",
                "6" = "#406792",
                "5" = "#B46F9C",
                "4" = "#D05B61",
                "3" = "#FE7B47",
                "2" = "#FFA72C",
                "1" = "#FFBD71")

#MoAs
MoAs_colors <- c("DNA related agent" = "#6cca8e", #This is the Pitiscale with some changes
                 "Cell cycle arrest" = "#8398dc",
                 "Combination drug" = "#ea95ae",
                 "PI3K/AKT/mTOR signaling inhibitor" = "#1dade6", 
                 "Microtubule agent" = "#ff5f76", 
                 "Chromatin agent" = "#ffb6b6",
                 "EGFR inhibitor" = "#fff154",
                 "Pro-apoptotic agent" = "#ba7fff",
                 "Kinase inhibitor" = "#ffdd72", 
                 "MAPK inhibitor" = "#4b71e5", 
                 "VEGFR inhibitor" = "#888888",
                 "Metabolism disruptor" = "#ff6600",
                 "HSP inhibitor" = "#add82f",
                 "Transcription inhibitor" = "#ff3333",
                 "IGF1R signaling inhibitor" = "#0dba3c", 
                 "Ubiquitin-proteasome system inhibitor" = "#ff864c", 
                 "p53 activator/MDM2 inhibitor" = "#c4ea94",
                 "ROS/RNS modulator" = "#666699",
                 "JAK-STAT signaling inhibitor" = "#454D66",
                 "BRAF inhibitor" = "#d58aca",
                 "Ferroptosis inducer" = "#6da753",
                 "NAMPT inhibitor" = "#ca9a8c",
                 "NFkB signaling inhibitor" = "#ff4430",
                 "ATP related agent" = "#e06d23",
                 "MET inhibitor" = "#b1ede8",
                 "BCR-ABL inhibitor" = "#8932a8",
                 "Protein synthesis inhibitor" = "#3267a8",
                 "Other" = "#BBB9B7")

## TME colors
tme_colors <- c(
    "Immune_stromal_cDC1_bias" = "#b6d7a8",
    "Immune_stromal_CAFlike_cDC1_cDC2_bias" = "#fce5cd",
    "Immune_rich_cDC1_bias" = "#c9daf8",
    "Immune_desert_cDC1_bias" = "#cdb4db",
    "Tcell_centric_Mo_cDC1_cDC2_bias" = "#ead1dc",
    "Immune_stromal_CAFlike" = "#fcf6bd",
    "Immune_desert_Endo_Mo_cDC2_bias" = "#dd7e6b",
    "Myeloid_centric" = "#d4a373",
    "Immune_stromal_Endolike_Treg_cDC1_cDC2_bias" = "#ffafcc",
    "Myeloid_centric_Mp_bias" = "#98f5e1",
    "Immune_rich_Treg_bias" = "#f08080",
    "Myeloid_centric_Mo_cDC1_cDC2_bias" = "#a2d2ff",
    "none" = "#BBB9B7"
    )

## Broad Cancer type colors.
broad_cancer_type_colors <- c("Blood Cancer" = "#A3181B",
                              "Brain Cancer" = "#B2509E",
                              "Neuroblastic Tumors" = "#FA8528",
                              "Breast Cancer" = "#db447a",
                              "Skin Cancer" =  "#5E2D2C",
                              "Lung Cancer" = "#158A88",
                              "Sarcoma/Soft Tissue Cancer" = "#e8d52c",
                              "Esophageal Cancer" = "#007EB5",
                              "Bladder Cancer" = "#367040",
                              "Liver/Biliary Cancer" = "#03543C",
                              "Pancreatic Cancer" = "#694E85",
                              "Ovarian Cancer" = "#E834EB",
                              "Prostate Cancer" = "#005D95",
                              "Colon/Colorectal Cancer" = "#a7495a",
                              "Endometrial/Uterine Cancer" = "#FAD2D9",
                              "Head and Neck Cancer" = "#97D1A9",
                              "Kidney Cancer" = "#918050",
                              "Other" = "#6b6363"
)
