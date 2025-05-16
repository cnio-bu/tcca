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
                      "Untreated" = "#D18B6E",
                      "Unknown" = "#BBB9B7")

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

# # MPs clusters x10
# mps_colors <-  c("MP1"  = "#FFBD71",
#                  "MP2"  = "#FFA72C",
#                  "MP3"  = "#FE7B47",
#                  "MP4"  = "#D05B61",
#                  "MP5" = "#FFB0EA",
#                  "MP6"  = "#B46F9C",
#                  "MP7" = "#A52390",
#                  "MP8" = "#6567BD",
#                  "MP9"  = "#406792",
#                  "MP10"  = "#369CBB",
#                  "MP11"  = "#90E2ED",
#                  "MP12"  = "#43978D",
#                  "MP13" = "#A7D676",
#                  "MP14" = "#D2E295",
#                  "MP15" = "#DBDF00",
#                  "MP16" = "#FFE364")

# scTherapy clusters
sctherapy_colors <- c("1" = "#FFBD71",
                      "2" = "#FFA72C",
                      "3" = "#FE7B47",
                      "4" = "#D05B61",
                      "5" = "#FFB0EA",
                      "6" = "#B46F9C",
                      "7" = "#A52390",
                      "8" = "#6567BD",
                      "9" = "#369CBB",
                      "10" = "#A7D676")

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
                 "Cell cycle arrest;PI3K/AKT/mTOR signaling inhibitor" = "#086AFC",
                 "SRC inhibitor" = "#FFFE89",
                 "Other" = "#BBB9B7")

## TME colors
tme_colors <- c("Immune_stromal_cDC1_bias" = "#006a40",
                "Immune_stromal_CAFlike_cDC1_cDC2_bias" = "#75b41e",
                "Immune_rich_cDC1_bias" = "#b3d0c0",
                "Immune_desert_cDC1_bias" = "#007e7f",
                "Tcell_centric_Mo_cDC1_cDC2_bias" = "#8ba1bc",
                "Immune_stromal_CAFlike" = "#f2990c",
                "Immune_desert_Endo_Mo_cDC2_bias" = "#d35e60",
                "Myeloid_centric" = "#f08892",
                "Immune_stromal_Endolike_Treg_cDC1_cDC2_bias" = "#8ab8cf",
                "Myeloid_centric_Mp_bias" = "#b4a814",
                "Immune_rich_Treg_bias" = "#9b7ebd",
                "Myeloid_centric_Mo_cDC1_cDC2_bias" = "#5a5895",
                "none" = "#BBB9B7")

## Broad Cancer type colors.
broad_cancer_type_colors <- c("Blood Cancer" = "#A3181B",
                              "Brain Cancer" = "#B2509E",
                              "Neuroblastic Tumors" = "#F06616",
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
                              "Miscellaneous Cancer" =  "#6b6363",
                              "Other" = "#6b6363"
)

# Study colors
study_colors <- c("adrenalnb_rui_chong" = "#75AE25",
                  "all_maxime_caron" = "#AFDB57",
                  "aml_audrey_lasry" = "#c3d680",
                  "aml_sander_lambo" = "#bbd7b9",
                  "bcc_catherine_dyao" = "#5f8772",
                  "bone_yun_liu" = "#8a7c31",
                  "brca_bhupinder_pal" = "#1679CD",
                  "breast_sunny_wu" = "#6978d0",
                  "brmets_hugo_gonzalez" = "#7bbdd2",
                  "brmets_jana_biermann" = "#5FE7F9",
                  "cc_xiaosong_lu" = "#628FFE",
                  "cell_lines_gabriella_kinker" = "#064D5B",
                  "chol_min_zhang" = "#7b3981",
                  "cll_ramon_massoni" = "#A96EB0",
                  "crc_florian_uhlitz" = "#9aa5da",
                  "eac_thomas_carroll" = "#322544",
                  "esca_xiannian_zhang" = "#d1a93c",
                  "gbm_nourhan_abdelfattah" = "#e38b2d",
                  "luad_kim_nayoung" = "#FDDA0D",
                  "luad_philip_bisschof" = "#CFD646",
                  "mmieloma_stephan_tirier" = "#FFEE8C",
                  "nsclc_stefan_salcher" = "#8B3100",
                  "oc_ec_matthew_regner" = "#9E5951",
                  "pancancer_dalia_barkley" = "#D30017",
                  "pancancer_junbin_qian" = "#de4330",
                  "pancancer_sunny_wu" = "#d63f76",
                  "pdac_junya_peng" = "#DA5C8C",
                  "pdac_shu_zhang" = "#d688c3",
                  "pleural_rui_dong" = "#950443",
                  "prad_sujun_chen" = "#da7772",
                  "rcell_kevin_bi" = "#c9a2ac",
                  "rcell_r_li" = "#8a3c5a",
                  "skcm_chao_zhang" = "#c16b38",
                  "synovial_jerby_arnon" = "#d6af81",
                  "urothelial_chen" = "#482323",
                  "uvm_michael_durante" = "#556785")

tumor_type_colors <-  c("ALL"   = "#1f77b4", 
                        "LAML"  = "#ff7f0e",
                        "BCC"   = "#2ca02c",  
                        "SARC"  = "#8B0000",
                        "BRCA"  = "#9467bd",
                        "SKCM"  = "#8c564b",
                        "SCLC"  = "#e377c2",
                        "NSCLC" = "#7f7f7f",
                        "OV"    = "#bcbd22",
                        "COAD"  = "#17becf",
                        "MISC"  = "#393b79",
                        "CESC"  = "#637939",
                        "CHOL"  = "#8c6d31",
                        "CLL"   = "#843c39",
                        "READ"  = "#7b4173",
                        "ESCA"  = "#1b9e77",
                        "GBM"   = "#d95f02",
                        "OGD"   = "#7570b3",
                        "LUAD"  = "#e7298a",
                        "MM"    = "#66a61e",
                        "LUSC"  = "#e6ab02", 
                        "LCLC"  = "#a6761d", 
                        "UCEC"  = "#666666", 
                        "GIST"  = "#1c9099", 
                        "PAAD"  = "#d73027",  
                        "LIHC"  = "#74add1", 
                        "PRAD"  = "#f46d43",  
                        "KIRC"  = "#a6d96a",  
                        "HNSC"  = "#fdae61",  
                        "SKSC"  = "#3288bd", 
                        "PLEU"  = "#abdda4",
                        "SKAM"  = "#5e4fa2",
                        "BLCA"  = "#9e0142",
                        "UVM"   = "#ffffbf" 
)
# ['#ffccff', '#000000', '#ffcc99', '#ff0099', '#330000',
#   '#ffcc00', '#00ffff', '#cccc00', '#ff9900', '#000066',
#   '#9999ff', '#cc0000', '#339999', '#cc99cc', '#3366ff',
#   '#cc9966', '#003399', '#009933', '#00ccff']
