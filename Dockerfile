# Use an official R base image
FROM r-base:4.3.1

# Install necessary system libraries
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff5-dev \
    pandoc \
    procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Change working directory
WORKDIR /futuriandgeDownstream

# Install your package from GitHub
RUN R -e "install.packages('remotes')"

# Install Bioconductor packages
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('rtracklayer', 'ReactomePA', 'org.Mm.eg.db', 'org.Hs.eg.db', 'EnhancedVolcano', 'DOSE', 'DESeq2', 'ComplexHeatmap', 'biomaRt'), ask = FALSE)"

# Install main package
RUN echo "Cache busting value: $(date)" && R -e "remotes::install_github('ajandria/futuriandgeDownstream', force = TRUE)"

# Copy test files into Dockerfile
COPY used-locally-for-testing used-locally-for-testing

RUN R -e "install.packages('tidyverse')"
RUN R -e "BiocManager::install('clusterProfiler', ask = FALSE)"
