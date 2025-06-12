# # Annotations table rendered and displayed:

# 1) UI:  
mod_annotations_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Dropdown to filter by study (allows multiple selections)
    selectInput(
      ns("study_filter"),
      label   = "Select study",
      choices = NULL,
      multiple = TRUE
    ),
    br(),
    # Interactive data table
    DTOutput(ns("annotations_table")),
    br(),
    # Download buttons
    downloadButton(ns("download_tsv"), "Download TSV"),
    downloadButton(ns("download_pdf"), "Download PDF")
  )
}

# 2) Server:
mod_annotations_server <- function(id, tsv_data, annotations_tsv_path) {
  moduleServer(id, function(input, output, session) {
    
    
    # 2.1) Populate study dropdown when data is available
    observeEvent(tsv_data(), {
      studies <- sort(unique(tsv_data()$study))
      updateSelectInput(
        session,
        "study_filter",
        choices  = studies,
        selected = studies
      )
    }, ignoreNULL = FALSE)
    
    # 2.2) Reactive filtered data
    filtered_data <- reactive({
      req(tsv_data(), input$study_filter)
      subset(tsv_data(), study %in% input$study_filter)
    })
    
    # 2.3) Render interactive data table
    output$annotations_table <- renderDT({
      datatable(
        filtered_data(),
        filter  = "top",
        options = list(pageLength = 10, scrollX = TRUE)
      )
    })
    
    # 2.4) Download handler for TSV file
    output$download_tsv <- downloadHandler(
      filename = function() {
        paste0("annotations_",
               paste(input$study_filter, collapse = "_"),
               ".tsv")
      },
      content = function(file) {
        write.table(
          filtered_data(),
          file,
          sep       = "\t",
          row.names = FALSE,
          quote     = FALSE
        )
      }
    )
    
    # 2.5) Download handler for PDF
    output$download_pdf <- downloadHandler(
      filename = function() {
        paste0("annotations_",
               paste(input$study_filter, collapse = "_"),
               ".pdf")
      },
      content = function(file) {
        # Open PDF device
        grDevices::pdf(file, width = 11, height = 8.5)
        
        # Draw table (or just head() if dataset is huge)
        gridExtra::grid.table(
          filtered_data(),
          rows = NULL
        )
        
        # Close the PDF device
        grDevices::dev.off()
      }
    )
    
  })
}