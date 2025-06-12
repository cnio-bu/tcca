# setup_shiny_app.R

# 0) Configurating Reticulate to use my environment
if (requireNamespace("reticulate", quietly = TRUE)) {
  reticulate::use_condaenv("r_shiny_env", required = TRUE)
} else {
  stop("El paquete reticulate no está instalado")
}

# 1) Install ShinyCell from GitHub if it's not alrasdy installed:
if (!requireNamespace("ShinyCell", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("SGDDNB/ShinyCell")
}

# 2) Install SeuratDisk from GitHub if it's not alrasdy installed:
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github("mojaveazure/seurat-disk")
}

# 3) Install 'anndata' Python package if it's not already installed
if (!"anndata" %in% reticulate::py_list_packages()$package) {
  reticulate::py_install("anndata", pip = TRUE)
}

cat("Dependencies setup: ShinyCell, SeuratDisk & anndata installed\n")