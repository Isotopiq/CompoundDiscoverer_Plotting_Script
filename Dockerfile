FROM rocker/shiny:4.5.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    ca-certificates \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c(\
  'shiny','readxl','dplyr','tidyr','stringr','ggplot2','ggrepel','pheatmap','DT',\
  'pROC','randomForest','e1071','pls','httr','jsonlite'\
), repos='https://cloud.r-project.org/')"

WORKDIR /srv/shiny-server
COPY app.R /srv/shiny-server/app.R

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
