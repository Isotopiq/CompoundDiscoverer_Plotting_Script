# Compound Discoverer Modern Shiny Bundle

This bundle contains a modern light dashboard wrapper around the same analysis engine.

Included:
- PCA
- volcano plot
- clustered heatmap
- sample summaries
- top differential table
- individual metabolite bar plots
- ZIP export
- progress bar during analysis

## Local R install
```r
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
```

## Docker / Easypanel
```bash
docker compose up --build -d
```

Then open:
http://YOUR_SERVER_IP:3838
