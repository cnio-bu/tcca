# CLINICAL METADATA module:
# We will keep the user interface and server functions in the same script for now):
mod_clinical_metadata_ui <- function(id) {
  ns <- NS(id)
  tagList(
    downloadButton(ns("dl"), "Download Clinical TSV"),
    br(), br(),
    DTOutput(ns("tbl"))
  )
}

mod_clinical_metadata_server <- function(id, data, path) {
  moduleServer(id, function(input, output, session) {
    output$tbl <- renderDT({
      datatable(data, filter="top", options=list(pageLength=10, scrollX=TRUE))
    })
    output$dl <- downloadHandler(
      filename = function() basename(path),
      content  = function(file) file.copy(path, file)
    )
  })
}