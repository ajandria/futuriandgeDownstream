# Use an official R base image
FROM r-base:4.3.1

# Install necessary system libraries
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff5-dev \
    pandoc \
    procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# First, you’ll need to get renv installed on your Docker image. The easiest way to accomplish this is with the remotes package. 
# For example, if you wanted to install a specific version of renv from GitHub:
RUN R -e "install.packages('remotes')"
RUN R -e "install.packages('renv')"

# Change working directory
WORKDIR /futuriandgeDownstream

# Install your package from GitHub
RUN R -e "install.packages('remotes')"
RUN echo "Cache busting value: $(date)" && R -e "remotes::install_github('ajandria/futuriandgeDownstream')"

# Copy test files into Dockerfile
COPY used-locally-for-testing used-locally-for-testing
