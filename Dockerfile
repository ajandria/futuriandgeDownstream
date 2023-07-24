# Use an official R base image
FROM r-base:4.3.1

# Install necessary system libraries
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev 

# Install your package's dependencies
RUN R -e "install.packages(c('readr', 'purrr', 'dplyr', 'testthat', 'devtools', 'roxygen2'))"