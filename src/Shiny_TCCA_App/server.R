# Upload file paths and pre-uploaded data in global.R:
#uploaded_files <- mod_upload_server("upload_module")

# Import modules for the UI intended functionalities:
source("modules/mod_clinical_metadata.R")
source("modules/mod_barplot_umap.R")
source("modules/mod_subclone.R")
source("modules/mod_downloads.R")
source("modules/mod_annotations.R")
source("global.R")


# Backend server logic: 
server <- function(input, output, session) {
  
  # Clinical Metadata module
  mod_clinical_metadata_server(
    "clinical",
    data = clinical_md,
    path = clinical_tsv_path
  )
  
  
  # Barplots and UMAPs linked to ShinyCell app module:
  mod_barplot_umap_server(
    "viz"
  )
  
  # Annotations module:
  mod_annotations_server(
    "ann",
    tsv_data           = reactive(annotations_raw),
    annotations_tsv_path = reactive("www/tcca_annotation_raw.tsv")
  )
  
  # Subclone module:
  mod_subclone_server(
    "subclone",
    default_data = subclone_tsv
  )
  
  
  
  # Global downloads module:
  mod_downloads_server(
    "download",
    p_h5ad = h5ad_raw_path,
    p_ann  = annotations_tsv_path,
    p_clin = clinical_tsv_path,
    p_sub  = subclone_tsv_path
  )
  
}
