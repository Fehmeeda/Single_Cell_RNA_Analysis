# Use an official Ubuntu as a parent image
FROM ubuntu:latest

# Add the CRAN repository to the software sources list
RUN echo "deb http://cran.stat.ucla.edu/bin/linux/ubuntu focal-cran40/" | tee -a /etc/apt/sources.list.d/r-project.list

# Install gnupg package
RUN apt-get update && apt-get install -y gnupg

# Add the repository authentication key
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Update package lists
RUN apt-get update

# Install R and software to compile R add-on packages
RUN apt-get install -y r-base r-base-dev

# Install additional packages needed for R script
RUN apt-get install -y \
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
    libx11-dev \
    libhdf5-dev

# Install Bioconductor and required R packages
RUN R -e "install.packages(c('BiocManager', 'devtools', 'R.utils','hdf5r'))"
RUN R -e "BiocManager::install(c('Seurat', 'ggplot2', 'tidyverse', 'gridExtra', 'SeuratObject', 'patchwork', 'celldex', 'openxlsx','SingleR', 'pheatmap', 'org.Hs.eg.db', 'readxl','viridis'))"

# Create a directory for scripts and datasets
RUN mkdir /data

# Copy scripts and datasets
COPY GSE233615-SingleCellRNAseq-script2.R /data/
COPY GSE118127-GSE233615_RAW-dataset/ /data/GSE118127-GSE233615_RAW-dataset/

# Set the working directory to /data
WORKDIR /data/

# Run R script
CMD ["Rscript", "GSE233615-SingleCellRNAseq-script2.R"]
