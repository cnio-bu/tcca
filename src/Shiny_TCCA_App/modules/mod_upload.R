# Upload documents module (we will keep the user interface and server functions in the same script for now):

# 1) UI:
mod_upload_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("tsv_file"), "Upload .tsv file:", accept = ".tsv"),
    fileInput(ns("h5ad_file"), "Upload .h5ad file:", accept = ".h5ad"),
    fileInput(ns("subclone_file"), "Upload subclone .tsv file:", accept = ".tsv"),
    downloadButton(ns("download_h5ad"), "Download .h5ad file")
  )
}


# 2) Server:
mod_upload_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Upload .tsv format:
     tsv_data <- reactive({
      req(input$tsv_file)
      read.delim(input$tsv_file$datapath, stringsAsFactors = FALSE)
    })
    
    
     # Provide h5ad path for download & UMAP
     h5ad_path <- reactive({
       req(input$h5ad_file)
       input$h5ad_file$datapath
     })
     
     output$download_h5ad <- downloadHandler(
       filename = function() {
         basename(input$h5ad_file$name)
       },
       content = function(file) {
         file.copy(h5ad_path(), file)
       }
     )
     
     
     # Subclone data
     subclone_data <- reactive({
       req(input$subclone_file)
       read.delim(input$subclone_file$datapath, stringsAsFactors = FALSE)
     })
     
     
     return(list(
       tsv_data = tsv_data,
       h5ad_path = h5ad_path,
       subclone_data = subclone_data
     ))
  })
}

