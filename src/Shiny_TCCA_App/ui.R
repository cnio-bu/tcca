# FRONTEND: User interface (UI) for the application

# Import modules for the UI intended functionalities:
source("modules/mod_clinical_metadata.R")
source("modules/mod_barplot_umap.R")
source("modules/mod_subclone.R")
source("modules/mod_downloads.R")
source("modules/mod_annotations.R")
source("global.R")



ui <- fluidPage(
  tags$head(
    tags$link(rel="stylesheet", type="text/css", href="custom.css")
  ),
  
  navbarPage(
    title = div(
      tags$img(src="logo_cnio.svg",    height="40px"),
      tags$span("The Therapeutic Cancer Cell Atlas", class="navbar-title"),
      tags$img(src="logo_tcca.png",     height="40px")
    ),
    id = "main_nav", windowTitle = "TCCA",
    
    tabPanel("Clinical Metadata",
             mod_clinical_metadata_ui("clinical")
    ),
    
    tabPanel("Annotations",
             mod_annotations_ui("ann")
    ),
    
    tabPanel("Barplots / UMAPs",
             mod_barplot_umap_ui("viz")
    ),
    
    tabPanel("Drug prediction for subclones",
             mod_subclone_ui("subclone")
    ),
    
    tabPanel("Downloads",
             mod_downloads_ui("download")
    )
  )
)