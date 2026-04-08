install.packages(c(
  "shiny",
  "readxl",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "DT"
), repos = "https://cloud.r-project.org/")

shiny::runApp("app.R")
