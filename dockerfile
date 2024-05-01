# Use an official Ubuntu as a parent image
FROM ubuntu:latest

# Set environment variables for R installation
ENV DEBIAN_FRONTEND noninteractive

# Update indices and install required packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    pkg-config \
    libharfbuzz-dev \
    libproj-dev \
    libcairo2-dev \
    libfribidi-dev \
    libjpeg-dev \
    libx11-dev

# Add the signing key (by Michael Rutter) for the CRAN repository
RUN curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | apt-key add -

# Add the R 4.0 repository from CRAN
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

# Update package manager again
RUN apt-get update

# Install the latest version of R
RUN apt-get install -y r-base r-base-dev

# Install Bioconductor and required R packages
RUN R -e "install.packages(c('BiocManager', 'devtools', 'R.utils'), repos='https://cloud.r-project.org/')"
RUN R -e "BiocManager::install(c('Seurat', 'ggplot2', 'tidyverse', 'gridExtra', 'SeuratObject', 'patchwork', 'celldex', 'openxlsx', 'SingleR', 'pheatmap', 'org.Hs.eg.db', 'readxl','viridis'))"

# Install ggtree from GitHub
RUN R -e "devtools::install_github('YuLab-SMU/ggtree')"

# Create a directory for scripts and datasets
RUN mkdir /data

# Copy scripts and datasets
COPY GSE233615-SingleCellRNAseq-script2.R /data/
COPY GSE118127-GSE233615_RAW-dataset/ /data/GSE118127-GSE233615_RAW-dataset/

# Set the working directory to /data
WORKDIR /data/

# Run R script
CMD ["Rscript", "GSE233615-SingleCellRNAseq-script2.R"]
