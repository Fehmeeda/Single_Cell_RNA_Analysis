# Use an official Ubuntu as a parent image
FROM ubuntu:latest

# Set environment variables for R installation
ENV DEBIAN_FRONTEND noninteractive

# Install required packages
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
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
    wget 



# Install Bioconductor
RUN R -e "install.packages('BiocManager')"

# Run R command to install devtools package
RUN R -e "install.packages(c('devtools','R.utils'))"

# Run R command to install ggtree from GitHub
RUN R -e "devtools::install_github('YuLab-SMU/ggtree')"

# Install R packages using BiocManager
RUN R -e "BiocManager::install(c('oligo', 'GenomicRanges', 'Biostrings', 'SummarizedExperiment', 'MatrixGenerics', 'DelayedArray', 'oligoClasses', 'Biobase', 'multiClust', 'limma', 'EnhancedVolcano', 'diffcoexp', 'enrichplot', 'pathview', 'org.Hs.eg.db', 'pheatmap', 'ggplot2', 'amap', 'ggrepel', 'openxlsx', 'readxl', 'ggridges','pd.hg.u133.plus.2'))"

# Install clusterProfiler from GitHub
RUN R -e "devtools::install_github('YuLab-SMU/clusterProfiler')"

# Create a directory for all the folders and scripts
RUN mkdir /data

# Copy scripts into the /data directory
COPY scripts/ /data/

# Copy required files and datasets into the /data directory
COPY GSE40595_RAW/ /data/GSE40595_RAW/
COPY Required_files/ /data/Required_files/

# Set the working directory to /data
WORKDIR /data/

# Run R scripts
CMD ["Rscript", "GSE40595_Microarray_data_analysis.R"]
