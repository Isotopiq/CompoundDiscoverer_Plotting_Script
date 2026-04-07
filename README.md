# Compound Discoverer Shiny App + Easypanel Deployment Bundle

## Included files
- `app.R` — the Shiny app
- `docker-compose.yml` — Easypanel/Docker Compose deployment
- `Dockerfile` — builds the Shiny container
- `.dockerignore`
- `install_dependencies.R` — local R dependency installation helper

## Local R install
Run this inside R:

```r
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

shiny::runApp("app.R")
```

Or source the helper file:

```r
source("install_dependencies.R")
```

## Docker / Easypanel
Put all files in one project directory and deploy with:

```bash
docker compose up --build -d
```

### Why there is no fixed host port
The compose file uses:

```yaml
expose:
  - "3838"
```

This is usually best for Easypanel because Easypanel handles public routing and avoids host-port conflicts.

### If you want Docker to assign a random host port
Edit `docker-compose.yml` and replace the `expose` block with:

```yaml
ports:
  - target: 3838
    published: "0"
    protocol: tcp
    mode: host
```

That tells Docker to choose an ephemeral host port.

## Manual local Docker build
```bash
docker compose build
docker compose up -d
```

## Notes
- The container already installs the required R packages during build.
- For Easypanel, upload/import the project and let Easypanel build from the included `Dockerfile`.
