# Barplot and UMAP module from th ShinyCell Github app: 
# The ShinyCell sub app is embedded.
# We will keep the user interface and server functions in the same script for now):

# 1) UI
mod_barplot_umap_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$div(
      id = ns("shinycell_container"),
      style = "
        position: absolute;
        top: 90px;      /* ajusta este valor a la altura real de tu navbar */
        bottom: 0;
        left: 0;
        right: 0;
        overflow: hidden;
        ",
      #Inyect the iframe from the server
      uiOutput(ns("shinycell_iframe"))
    )
  )
}

# 2) SERVER 
mod_barplot_umap_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    output$shinycell_iframe <- renderUI({
      protocol <- session$clientData$url_protocol       # "http:" o "https:"
      hostname <- session$clientData$url_hostname       # "mi-servidor.tld" o "localhost"
      portMain <- session$clientData$url_port           # "3838" (app principal)
      # construimos la URL para puerto 5432
      url_subapp <- paste0(protocol, "//", hostname, ":19002")
      
      tags$iframe(
        src         = url_subapp,
        style       = "
          position: absolute;
          top: 0; bottom: 0;
          left: 0; right: 0;
          width: 100%;
          height: 100%;
          border: none;
        "
      )
      
    })
    
  })
} 