## TCCA COLOR PALETTES
# Cancer type colors 
# Tumor site x18 (17 + Other)
tumor_sites_colors <- c("Bone marrow" = "#FFE072",
                        "Brain" = "#CA9A8C",
                        "Adrenal gland" = "#FF9A4A",
                        "Breast" = "#FF8FAB",
                        "Skin" = "#5E2D2C",
                        "Lung" = "#A2F0EF",
                        "Soft tissue" = "#FB6467",
                        "Esophagus" = "#5780FE", 
                        "Lymph node" = "#B6F884",
                        "Liver" = "#309975",
                        "Pancreas" = "#B47EB3",
                        "Ovary" = "#9C0D38",
                        "Prostate" = "#005D95",
                        "Colon" = "#BD081A",
                        "Kidney"  = "#918050",
                        "Bladder" = "#694E85",
                        "Other" = "#BBB9B7")

# Cancer type.
cancer_type_colors <- list(
  "Brain Cancer" = "#CA9A8C",
  "Neuroblastic Tumors" = "#FF9A4A",
  "Blood Cancer" = "#FFE072",
  "Skin Cancer" = "#5E2D2C",
  "Sarcoma/Soft Tissue Cancer" = "#FB6467",
  "Breast Cancer" = "#FF8FAB",
  "Lung Cancer" = "#A2F0EF",
  "Ovarian Cancer" = "#9C0D38",
  "Colon/Colorectal Cancer" = "#BD081A",
  "Endometrial/Uterine Cancer" = "#F7C0AC",
  "Liver/Biliary Cancer" = "#309975",
  "Bladder Cancer" = "#694E85",
  "Head and Neck Cancer" = "#83343D",
  "Prostate Cancer" = "#005D95",
  "Kidney Cancer" = "#918050",
  "Esophageal Cancer" = "#5780FE",
  "Pancreatic Cancer" = "#B47EB3",
  "Thyroid Cancer" = "#01ADB9",
  "Gastric Cancer" = "#A8A5D1",
  "Miscellaneous Cancer" = "#BBB9B7")

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

# Cell line/patient
ccl_p_colors <- c("Patient" = "#427394", 
                 "Cell line" = "#F78C1F")

# Study
study_colors <- c(
  "adrenalnb_rui_chong" = "#304428",
  "all_maxime_caron" = "#682EDE",
  "aml_audrey_lasry" = "#63D945",
  "aml_sander_lambo" = "#CE38E1",
  "bcc_catherine_dyao" = "#A9DA41",
  "bone_yun_liu" = "#5F41C0",
  "brca_bhupinder_pal" = "#D7D23D",
  "breast_sunny_wu" = "#B668DD",
  "brmets_hugo_gonzalez" = "#65D888",
  "brmets_jana_biermann" = "#CC41AD",
  "cc_xiaosong_lu" = "#538A3D",
  "cell_lines_gabriella_kinker" = "#382776",
  "chol_min_zhang" = "#D1A93C",
  "cll_ramon_massoni" = "#6978D0",
  "crc_florian_uhlitz" = "#E38B2D",
  "eac_thomas_carroll" = "#7B3981",
  "esca_xiannian_zhang" = "#C3D680",
  "gbm_nourhan_abdelfattah" = "#D63F76",
  "luad_kim_nayoung" = "#60D6C0",
  "luad_philip_bisschof" = "#DE4330",
  "mmieloma_stephan_tirier" = "#7BBDD2",
  "nsclc_stefan_salcher" = "#8E3026",
  "oc_ec_matthew_regner" = "#BBD7B9",
  "pancancer_dalia_barkley" = "#322544",
  "pancancer_junbin_qian" = "#D6AF81",
  "pancancer_sunny_wu" = "#482323",
  "pdac_junya_peng" = "#9AA5DA",
  "pdac_shu_zhang" = "#C16B38",
  "pleural_rui_dong" = "#556785",
  "prad_sujun_chen" = "#8A7C31",
  "rcell_kevin_bi" = "#D688C3",
  "rcell_r_li" = "#5F8772",
  "skcm_chao_zhang" = "#DA7772",
  "synovial_jerby_arno" = "#7B5839",
  "urothelial_chen" = "#C9A2AC",
  "uvm_michael_durante" = "#8A3C5A"
)


# Therapeutic clusters x5
tcs_colors <- c("4"= "#F8766D",
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
                 "Other" = "#BBB9B7")


## Therapeutic module colors
module_colors <- c(
    "TM1" = "#b6d7a8",
    "TM2" = "#fce5cd",
    "TM3" = "#c9daf8",
    "TM4" = "#d9d2e9",
    "TM5" = "#ead1dc",
    )

## TME archeatype
tme_colors <- list(
    "Immune_stromal_cDC1_bias" = "#006A40",
    "none" = "#95828D",
    "Immune_stromal_CAFlike_cDC1_cDC2_bias" = "#75B41E",
    "Myeloid_centric" = "#F08892",
    "Immune_rich_cDC1_bias" = "#708C98",
    "Immune_stromal_Endolike_Treg_cDC1_cDC2_bias" = "#8AB8CF",
    "Immune_desert_cDC1_bias" = "#007E7F",
    "Myeloid_centric_Mp_bias" = "#358359",
    "Tcell_centric_Mo_cDC1_cDC2_bias" = "#8BA1BC",
    "Immune_rich_Treg_bias" = "#9B7EBD",
    "Immune_stromal_CAFlike" = "#F2990C",
    "Myeloid_centric_Mo_cDC1_cDC2_bias" = "#5A5895",
    "Immune_desert_Endo_Mo_cDC2_bias" = "#D35E60"
)