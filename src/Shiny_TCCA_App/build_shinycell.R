## Building the github app ShinyCell:
library(ShinyCell)


# Path to .h5ad file (which will be our data source to print the UMAPS and barplots)
inpFile <- "www/tcca_annotated_clustered.h5ad"

# Shiny cell app configuration
scConf <- createConfig(inpFile)

# Shiny Cell app generation inside /www/shinycell_app
makeShinyApp(
  obj     = inpFile,
  scConf      = scConf,
  shiny.dir   = "shinyAppH5ad",
  shiny.title = "TCCA Single‑Cell Viewer"
)

