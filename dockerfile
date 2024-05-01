# Use an official Ubuntu as a parent image
FROM ubuntu:latest

# Set environment variables for R installation
ENV DEBIAN_FRONTEND noninteractive

# Update indices
RUN apt update -qq

# Install two helper packages we need
RUN apt install --no-install-recommends -y software-properties-common dirmngr

# Install curl
RUN apt-get update && apt-get install -y curl

# Add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

# Add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

# Update package manager again
RUN apt update

# Install the new version of R
RUN apt install -y r-base r-base-dev

# Install required packages
RUN apt-get install -y \ 
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    pkg-config\
    libharfbuzz-dev \
    libproj-dev \
    libcairo2-dev \
    libfribidi-dev \
    libjpeg-dev \
    libx11-dev 

# Install Bioconductor
RUN R -e "install.packages('BiocManager')"

# Run R command to install devtools package
RUN R -e "install.packages(c('devtools','R.utils'))"

# Run R command to install ggtree from GitHub
RUN R -e "devtools::install_github('YuLab-SMU/ggtree')"

# Install R packages using BiocManager
RUN R -e "BiocManager::install(c('Seurat', 'ggplot2', 'tidyverse', 'gridExtra', 'SeuratObject', 'patchwork', 'celldex', 'openxlsx', 'SingleR', 'pheatmap', 'org.Hs.eg.db', 'readxl','viridis'))"

# Create a directory for all the folders and scripts
RUN mkdir /data

# Copy scripts into the /data directory
COPY GSE233615-SingleCellRNAseq-script2.R /data/

# Copy required files and datasets into the /data directory
COPY GSE118127-GSE233615_RAW-dataset/ /data/GSE118127-GSE233615_RAW-dataset/

# Set the working directory to /data
WORKDIR /data/

# Run R scripts
CMD ["Rscript", "GSE233615-SingleCellRNAseq-script2.R"]
