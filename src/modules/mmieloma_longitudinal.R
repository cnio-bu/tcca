library(igraph)
library(tidyverse)

mmieloma_communities <- data.table::fread(
    "results/modules/annotated/MM_treated_communities.tsv"
    )

clinical_data = data.table::fread("results/annotation/clinical_metadata_v3_clean.tsv")


mielomas <- clinical_data %>%
    filter(tumor_type == "MM", treated == "t" )

rm006 <- clinical_data %>%
    filter(patient == "RRMM06")

samples_to_use <- rm006 %>%
    pull(sample)

drug_data <- data.table::fread("reference/final_moas - Collapsed.tsv")

## WHICH COMMS ARE RRMM06
rm006_communities <- mmieloma_communities %>%
    filter(sample %in% samples_to_use) %>%
    group_by(sample) %>%
    reframe(communities = community) %>%
    distinct()


comm_summary <- mmieloma_communities %>%
    group_by(community, signature) %>%
    mutate(
        recurrency = n()
    ) %>%
    filter(community == 5) %>%
    left_join(y = drug_data[,c("IDs", "studies", "preferred.drug.names")],
              by = c("signature" = "IDs"), multiple = "first"
              )


## Community 5 top Hit is a ROSN inhibitor, known to happen in MM (ROS)
## It is also prominent in RM006 T3. Habria que ver la historia BRD-K71935468


#### NECESITO
###
### Tratamientos anotados del paciente
### Lineas de respuesta
### beyondcell completo para rrmm006 aunque NO esten todas las celulas en bc-meta

#Preguntas: 
#    1. Que modulos estan representados en cada sample de RRMM066
#    2. Miramos las drogas. Casan ?
#    3. Superimponer los modulos.
#    4. Superimponer metacomunidades



