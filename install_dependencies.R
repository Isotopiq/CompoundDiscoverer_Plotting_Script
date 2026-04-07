# Install required R packages
install.packages(c(
  "shiny",
  "shinythemes",
  "readxl",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "DT"
), repos = "https://cloud.r-project.org/")

# Run locally
shiny::runApp("app.R")
