# Subclone module:
# We will keep the user interface and server functions in the same script for now):

# 1) UI
mod_subclone_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Raw Subclone Data:"),
    DTOutput(ns("subclone_table")),
    downloadButton(ns("download_full"), "Download full TSV"),
    hr(),
    h4("Drugs per Cancer Type:"),
    selectInput(ns("cancer_type"), "Select cancer type:", choices = NULL),
    DTOutput(ns("drugs_by_cancer")),
    downloadButton(ns("download_grouped"), "Download grouped TSV")
  )
}



# 2) SERVER:
mod_subclone_server <- function(id, default_data) {
  moduleServer(id, function(input, output, session) {
    
    # Update cancer type choices
    observe({
      types <- sort(unique(default_data$Broad.Tumor.Type))
      updateSelectInput(session, "cancer_type", choices = types, selected = types[1])
    })
    
    # Render full subclone table
    output$subclone_table <- renderDT({
      datatable(default_data, filter = "top", options = list(scrollX = TRUE, pageLength = 10))
    })
    
    # Download full table
    output$download_full <- downloadHandler(
      filename = function() "subclone_full.tsv",
      content = function(file) {
        write.table(default_data, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    
    # Grouped drug list by selected cancer type
    grouped <- reactive({
      req(input$cancer_type)
      default_data %>%
        filter(Broad.Tumor.Type == input$cancer_type) %>%
        distinct(Drug.Name) %>%
        arrange(Drug.Name)
    })
    
    # Render grouped table
    output$drugs_by_cancer <- renderDT({
      datatable(grouped(), options = list(dom = 't', pageLength = 10))
    })
    
    # Download grouped table
    output$download_grouped <- downloadHandler(
      filename = function() paste0("drugs_for_", input$cancer_type, ".tsv"),
      content = function(file) {
        write.table(grouped(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
  })
}