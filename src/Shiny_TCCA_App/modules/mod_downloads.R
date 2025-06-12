# Donwload files module: 
# We will keep the user interface and server functions in the same script for now):

# 1) UI:
mod_downloads_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Download Data Files"),
    downloadButton(ns("dl_h5ad"),    "Download raw .h5ad"),
    downloadButton(ns("dl_ann"),     "Download annotations .tsv"),
    downloadButton(ns("dl_clin"),    "Download clinical .tsv"),
    downloadButton(ns("dl_sub"),     "Download subclone .tsv")
  )
}



# 2) Server:
mod_downloads_server <- function(id, p_h5ad, p_ann, p_clin, p_sub) {
  moduleServer(id, function(input, output, session) {
    output$dl_h5ad <- downloadHandler(
      filename = function() basename(p_h5ad),
      content  = function(file) file.copy(p_h5ad, file)
    )
    output$dl_ann <- downloadHandler(
      filename = function() basename(p_ann),
      content  = function(file) file.copy(p_ann, file)
    )
    output$dl_clin <- downloadHandler(
      filename = function() basename(p_clin),
      content  = function(file) file.copy(p_clin, file)
    )
    output$dl_sub <- downloadHandler(
      filename = function() basename(p_sub),
      content  = function(file) file.copy(p_sub, file)
    )
  })
}