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
    libtiff5-dev

# First, you’ll need to get renv installed on your Docker image. The easiest way to accomplish this is with the remotes package. 
# For example, if you wanted to install a specific version of renv from GitHub:
RUN R -e "install.packages('remotes')"
RUN R -e "install.packages('renv')"

# Change working directory
WORKDIR /futuriandgeDownstream
COPY renv.lock renv.lock

# Next, you need to tell renv which library paths to use for package installation. 
# You can either set the RENV_PATHS_LIBRARY environment variable to a writable path within your Docker container, 
# or copy the renv auto-loader tools into the container so that a project-local library can be automatically provisioned and used when R is launched.
# approach two
RUN mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Finally, you can run renv::restore() to restore packages as defined in the lockfile:
RUN R -e "renv::restore()"

# Install your package from GitHub
RUN R -e "install.packages('remotes')"
RUN R -e "remotes::install_github('ajandria/futuriandgeDownstream', ref = 'HEAD')"

# Copy test files into Dockerfile
COPY used-locally-for-testing used-locally-for-testing
