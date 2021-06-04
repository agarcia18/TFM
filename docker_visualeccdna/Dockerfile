# Dockerfile for visualeccdna

FROM rocker/shiny:4.1.0

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libxml2 \
    libxml2-dev \
    r-cran-rstan \
    r-cran-xml
  
## rm -rf /var/lib/apt/lists/*
# install R packages required 
RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('rstan', repos='http://cran.rstudio.com/')"
RUN R -E "install.packages('igraph',repos='http://cran.rstudio.com/')"

# bioC 3.13
RUN R -e "BiocManager::install(version = '3.13', ask = FALSE)"
RUN R -e "devtools::install_version('shiny', version = '1.6.0', repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_version('shinydashboard', version = '0.7.1', repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_version('dplyr', version = '1.0.6', repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_version('ggplot2', version = '3.3.3', repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_version('stringr', version = '1.4.0', repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_version('plotly', version = '4.9.3', repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_version('DT', version = '0.18', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('TxDb.Hsapiens.UCSC.hg38.knownGene', 'GenomicRanges', 'org.Hs.eg.db', 'RITANdata', 'Biostrings','BSgenome.Hsapiens.UCSC.hg38'), version = '3.13')"
RUN R -e "BiocManager::install(c('RITAN'), version = '3.13')"
# code

COPY visualeccdna ./app

# COPY R /srv/shiny-server/R

# select port
EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]