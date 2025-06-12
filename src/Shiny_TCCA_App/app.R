# 1) Paths to the ShinyCell subapp and main app:
shinycell_subapp_path <- "/Users/isanzp/Documents/CNIO/Shiny/Shiny_TCCA_App/shinyAppH5ad"
main_app_path <- "/Users/isanzp/Documents/CNIO/Shiny/Shiny_TCCA_App"



# 2) Rebuild ShinyCell app if it doesn't exist
if (!dir.exists(shinycell_subapp_path) ||
    !file.exists(file.path(shinycell_subapp_path, "ui.R")) ||
    !file.exists(file.path(shinycell_subapp_path, "server.R"))) {
  source("build_shinycell.R")
}


# 3) Launch ShinyCell subapp in the background:
message("→ Launching ShinyCell sub-app on port 19002…")
cmd <- sprintf(
  "shiny::runApp('%s', port = 19002, host = '0.0.0.0', launch.browser = FALSE)",
  shinycell_subapp_path
)
processx::process$new(
  "Rscript",
  c("-e", cmd),
  stdout = "|",
  stderr = "2>&1",
  supervise = TRUE
)

# 4) Launch the main app on the specific port (19003):
message("→ Launching main app on port 19003")
shiny::runApp(
  appDir         = ".",
  host           = "0.0.0.0",
  port           = 19003,
  launch.browser = TRUE
)