## TCCA COLOR PALETTES

# Tumor type x15 (14 + Other)
tumor_sites_colors <- c(
  "Bone marrow" = "#FB6467",
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
  "Other" = "#6b6363"
)

# Patient/cell line
patient_ccl_colors <- c("Patient" = "#427394", "Cell line" = "#F78C1F")

# Treated/Untreated
treatment_colors <- c(
  "Treated" = "#6ED1BC",
  "Untreated" = "#D18B6E",
  "Unknown" = "#BBB9B7"
)

# Adult/Pediatric
age_colors <- c("Adult" = "#7689DE", "Pediatric" = "#a9dce3")
age_group_colors <- c(
  "Pediatric"   = "#a9dce3",
  "Young adult" = "#6CC2BD",
  "Adult"       = "#7689DE",
  "Elderly"     = "#4B4A73",
  "Unknown"     = "#BEBEBE"
)

# Male/Female
sex_colors <- c("Male" = "#ec5c44", "Female" = "#1c8c8c")

# Solid/Liquid
sl_colors <- c("Solid" = "#F6Bd60", "Liquid" = "#706695")

# Primary/Met
pm_colors <- c("Metastasis" = "#C10044",
               "Primary" = "#F0BFD0")

# Therapeutic clusters x5
tcs_colors <- c(
  "5" = "#A3A500",
  "4" = "#F8766D",
  "3" = "#FFD64C",
  "2" = "#00BF7D",
  "1" = "#00B0F6",
  "0" = "#E76BF3"
)

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
sctherapy_colors <- c(
  "1" = "#FFBD71",
  "2" = "#FFA72C",
  "3" = "#FE7B47",
  "4" = "#D05B61",
  "5" = "#FFB0EA",
  "6" = "#B46F9C",
  "7" = "#A52390",
  "8" = "#6567BD",
  "9" = "#369CBB",
  "10" = "#A7D676"
)

#MoAs
MoAs_colors <- c(
  "DNA related agent" = "#6cca8e",
  #This is the Pitiscale with some changes
  "Cell cycle arrest" = "#8398dc",
  "Combination drug" = "#ea95ae",
  "PI3K/AKT/mTOR signaling inhibitor" = "#1dade6",
  "Microtubule agent" = "#ff5f76",
  "Chromatin agent" = "#ffb6b6",
  "EGFR inhibitor" = "#fff154",
  "EGFR inhibitor;VEGFR inhibitor" = "#3f685e",
  "PARP inhibitor" = "#7d58ad",
  "Pro-apoptotic agent" = "#ba7fff",
  "Kinase inhibitor" = "#ffdd72",
  "MAPK inhibitor" = "#4b71e5",
  "VEGFR inhibitor" = "#6e4b3a",
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
  "Other" = "#BBB9B7"
)

## TME colors
tme_colors <- c(
  "Immune_stromal" = "#006a40",
  "Immune_stromal_CAFlike_Mp_bias" = "#75b41e",
  "Immune_rich" = "#b3d0c0",
  "Immune_stromal_desert" = "#007e7f",
  "Tcell_centric" = "#8ba1bc",
  "Immune_desert_CAFlike" = "#f2990c",
  "Myeloid_centric_Mo_bias(medium)" = "#d35e60",
  "Myeloid_centric" = "#f08892",
  "Immune_stromal_Endolike" = "#8ab8cf",
  "Myeloid_centric_Mp_bias" = "#b4a814",
  "Immune_rich_Treg_cDC2_bias" = "#9b7ebd",
  "Myeloid_centric_Mo_bias(high)" = "#5a5895",
  "none" = "#BBB9B7"
)


tme_group_colors <- c(
  "Immune_stromal" = "#1b9e77",
  "Immune_rich" = "#d95f02",
  "Immune_desert" = "#7570b3",
  "Tcell_centric" = "#e7298a",
  "Myeloid_centric" = "#41a24b",
  "none" = "#BBB9B7"
)

## Broad Cancer type colors.
broad_cancer_type_colors <- c(
  "Blood Cancer" = "#A3181B",
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
  "Gastric Cancer" = "#c3d9ff",
  "Thyroid Cancer" = "#ada631",
  "Other" = "#6b6363"
)

# Cell type colors
cell_type_colors <- c(
  "Stromal cell" = "#d7bafd",
  "Glial cell" = "#873e23",
  "Epithelial" = "#006cd1",
  "B-cell" = "#935ee7",
  "Dendritic cell" = "#009f83",
  "Plasmacytoid dendritic cell" = "#b3d266",
  "Endothelial" = "#d20086",
  "Erythrocyte" = "#7dda84",
  "Granulocyte" = "#ff44d4",
  "Mast" = "#0092f4",
  "Monocyte/Macrophage" = "#b66c00",
  "NK cell" = "#b88c9c",
  "Neuron" = "#cb4dd0",
  "Stem" = "#5fcaff",
  "Innate lymphoid cells" = "#7591ff",
  "Plasma cell" = "#f97d1b",
  "CD4+ T-cell" = "#394f94",
  "CD8+ T-cell" = "#d1da32",
  "Regulatory T-cell" = "#8e2c68",
  "Unconventional T-cells" = "#3d6700",
  "Malignant" = "#d82500",
  "Unknown" = "#808080"
)
# Study colors
study_colors <- c(
  "adrenalnb_rui_chong" = "#75AE25",
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
  "uvm_michael_durante" = "#556785"
)

tumor_type_colors <- c(
  "ALL"  = "#FF6600",
  "BCC"  = "#D62727",
  "BLCA" = "#FAD2D9",
  "BRCA" = "#ED2891",
  "CESC" = "#F6B667",
  "CHOL" = "#104A7F",
  "CLL"  = "#BD5D09",
  "COAD" = "#9EDDF9",
  "ESCA" = "#007EB5",
  "GBM"  = "#B2509E",
  "HNSC" = "#97D1A9",
  "KIRC" = "#F8AFB3",
  "LAML" = "#754C29",
  "LCLC" = "#98909E",
  "LGG"  = "#D49DC7",
  "LIHC" = "#CACCDB",
  "LUAD" = "#D3C3E0",
  "LUSC" = "#A084BD",
  "MESO" = "#542C88",
  "MM"   = "#BF4362",
  "NB"   = "#FFD700",
  "OV"   = "#D97D25",
  "PAAD" = "#6E7BA2",
  "PRAD" = "#7E1918",
  "READ" = "#DAF1FC",
  "SARC" = "#00A99D",
  "SCLC" = "#BAB095",
  "SKCM" = "#BBD642",
  "SKSC" = "#FF6666",
  "STAD" = "#00AEEF",
  "THCA" = "#F9ED32",
  "UCEC" = "#FBE3C7",
  "UCS"  = "#F89420",
  "UVM"  = "#009444"
)

mp_colors <- c(
  "MP1_CellCycle.G2M" = "#1b4248",
  "MP2_CellCycle.G1S" = "#559bab",
  "MP3_CellCycle.HMG-rich" = "#ffb0ea",
  "MP4_CellCycle.ChromatinRemodeling" = "#fe7b47",
  "MP5_CellCycle..Chromatin" = "#ffa72c",
  "MP6_CellCycle.DNARepair" = "#d05b61",
  "MP7_Oncogenic.MYC" = "#6567bd",
  "MP8_Stress.IEGs.AP1" = "#b46f9c",
  "MP9_Stress.ISR" = "#ffbd71",
  "MP10_Stress.Hypoxia" = "#a9b366",
  "MP11_Stress.Metabolic" = "#459d6b",
  "MP12_Stress.Detoxification" = "#00b6c8",
  "MP13_Stress.OxidativeStress" = "#953f6d",
  "MP14_Inflammation.Interferon.MHCII" = "#584116",
  "MP15_Inflammation.ReactiveEpithelia" = "#ea1755",
  "MP16_Inflammation.TNFA.NFkB" = "#8d0048",
  "MP17_EMT.partialEMT" = "#e4028c",
  "MP18_EMT.EMT_I" = "#ca9a8c",
  "MP19_EMT.EMT_II" = "#a554a1",
  "MP20_EMT.Mesenchymal-like" = "#ff8fab",
  "MP21_EMT.EMT_III" = "#bfdacb",
  "MP22_CellularPlasticity.ActiveSignaling" = "#c8d45e",
  "MP23_CellularPlasticity.Post-transcriptionalChromatin" = "#ffe700",
  "MP24_CellularPlasticity.EpithelialRemodeling" = "#f9b9b9",
  "MP25_ProteinRegulation.ProteasomeDegradation" = "#fcb451",
  "MP26_ProteinRegulation.ProteinMaturation" = "#937656",
  "MP27_ProteinRegulation.UPR" = "#8baaa6",
  "MP28_ProteinRegulation.ProteinTranslation" = "#5e91ab",
  "MP29_EpithelialSenescence" = "#6d8e43",
  "MP30_MitochondrialRespiration" = "#9d7e02",
  "MP31_Cilia" = "#0074bc",
  "MP32_LineageSpecific.Hemato.InflammatoryMyeloid" = "#4e5771",
  "MP33_LineageSpecific.Hemato.HPSCs" = "#ad8089",
  "MP34_LineageSpecific.Hemato.Neutrophil" = "#f28bb9",
  "MP35_LineageSpecific.Hemato.APC-MHCII" = "#f7951f",
  "MP36_LineageSpecific.Hemato.MastCells" = "#ccde66",
  "MP37_LineageSpecific.Hemato.RBCs" = "#3b803d",
  "MP38_LineageSpecific.Neural.Astrocytes" = "#95d4aa",
  "MP39_LineageSpecific.Neural.OPCs-NPCs" = "#be7d6d",
  "MP40_LineageSpecific.Melanocyte-Pigmentation" = "#56728a",
  "MP41_LineageSpecific.Prostate-Secretory" = "#93114f",
  "MP42_LineageSpecific.Urothelial-Secretory" = "#bd58a2",
  "MP43_LineageSpecific.Gastrointestinal-Secretory" = "#c07c2c"
)


mp_family_colors <- c(
  "CellCycle" = "#1b4248", 
  "Oncogenic" = "#6567bd",
  "Stress" = "#b46f9c",
  "Inflammation" = "#584116",
  "EMT" = "#e4028c", 
  "CellularPlasticity" = "#c8d45e",
  "ProteinRegulation" = "#fcb451",
  "EpithelialSenescence" = "#6d8e43",
  "MitochondrialRespiration" = "#9d7e02",
  "Cilia" = "#0074bc",
  "LineageSpecific.Hemato" = "#4e5771",
  "LineageSpecific.Neural" = "#95d4aa", 
  "LineageSpecific.Other" = "#9b2727ff"
)
