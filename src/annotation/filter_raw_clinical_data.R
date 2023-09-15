library("tidyverse")

raw_clinical_export <- data.table::fread("results/annotation/clinical_metadata_v1.tsv")

## Remove duplicated samples
samples_to_remove <- c(
"Kim_Lee_2020_BRONCHO_11",
"Kim_Lee_2020_BRONCHO_58",
"Kim_Lee_2020_EBUS_06",
"Kim_Lee_2020_EBUS_10",
"Kim_Lee_2020_EBUS_13",
"Kim_Lee_2020_EBUS_19",
"Kim_Lee_2020_EBUS_28",
"Kim_Lee_2020_EBUS_51",
"Kim_Lee_2020_LUNG_T18",
"Kim_Lee_2020_LUNG_T20",
"Kim_Lee_2020_LUNG_T28",
"Kim_Lee_2020_LUNG_T34",
"Kim_Lee_2020_NS_03",
"Kim_Lee_2020_NS_04",
"Kim_Lee_2020_NS_07",
"Kim_Lee_2020_NS_12",
"Kim_Lee_2020_NS_17",
"Kim_Lee_2020_NS_19",
"CHCAHsZhang18",
"CHCAHsZhang20",
"CHCAHsZhang23",
"CHCAHsZhang241",
"CHCAHsZhang242",
"LUADHsKimT18",
"LUADHsKimT20",
"LUADHsKimT30",
"LUADHsKimT34",
"BT1290",
"BT1291",
"BT1292",
"BT1295",
"BT1296",
"BT1297",
"BT1298",
"BT1299",
"BT1300"
)

clinical_unique_samples <- raw_clinical_export %>%
    filter(
        !(sample %in% samples_to_remove) 
    ) %>%
    mutate(
        refined_tumor_site = case_when(
            tumor_site == "adrenal gland" ~ "adrenal_gland",
            tumor_site == "cecum" ~ "colon",
            tumor_site == "colon adjacent" ~ "colon",
            tumor_site == "Face" ~ "skin",
            tumor_site == "Finger" ~ "skin",
            tumor_site == "haematopoietic_and_lymphoid_tissue" ~ "lymph_node",
            tumor_site == "Head" ~ "skin",
            tumor_site == "kidney adjacent" ~ "kidney",
            tumor_site == "knee" ~ "skin",
            tumor_site == "lung adjacent" ~ "lung",
            tumor_site == "mucosa adjacent" ~ "bladder",
            tumor_site == "Neck" ~ "skin",
            tumor_site == "pancreas adjacent" ~ "pancreas",
            tumor_site == "para_aortic" ~ "soft_tissue",
            tumor_site == "rectum_adjacent" ~ "rectum",
            tumor_site == "Sole" ~ "skin",
            tumor_site == "subcutaneous" ~ "soft_tissue",
            tumor_site == "thigh" ~ "skin",
            tumor_site == "tongue" ~ "oral_cavity",
            tumor_site == "uriry_tract" ~ "urinary_tract",
            tumor_site == "abdomen" ~ "abdomen_wall",
            TRUE ~ tumor_site
            
        )
    )

write_tsv(
    x = clinical_unique_samples,
    file = "results/annotation/clinical_metadata_v1_clean.tsv"
    )
