FROM jupyter/datascience-notebook

# ####################################
# INSTALLING R & packages
# ####################################
USER root

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    r-base \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \ 
    lmodern \
    default-jre

# Install required R packages
# RUN R -e "install.packages('Seurat', repos='https://cran.r-project.org/')"  
# RUN R -e "install.packages('remotes', repos='https://cran.r-project.org/')" 
# RUN R -e "install.packages('anndata', repos='https://cran.r-project.org/')" 
# RUN R -e "install.packages('devtools',repos='https://cran.r-project.org/')"
# RUN R -e "devtools::install_github('satijalab/seurat-data')"
# RUN R -e "install.packages('BiocManager', repos='https://cran.r-project.org/')" 
# RUN R -e "BiocManager::install('SingleCellExperiment')" 
# RUN R -e "BiocManager::install('scry')" 
# RUN R -e "BiocManager::install('DESeq2')"

WORKDIR /APP

# Copy requirements and install packages
#RUN pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
COPY ./requirements.txt /APP
RUN pip install --no-cache-dir -r requirements.txt

# 
# Declare web port used by jupyter lab

# RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.3.450/quarto-1.3.450-linux-amd64.deb && \
#     dpkg -i quarto-1.3.450-linux-amd64.deb && \
#     rm quarto-1.3.450-linux-amd64.deb

RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.7.33/quarto-1.7.33-linux-amd64.deb && \
    dpkg -i quarto-1.7.33-linux-amd64.deb && \
    rm quarto-1.7.33-linux-amd64.deb


RUN mkdir /tmp/numba_cache \
    && chmod 777 /tmp/numba_cache
ENV NUMBA_CACHE_DIR /tmp/numba_cache

EXPOSE 8989

# Set default command
#CMD start.sh jupyter lab --port=8889 --no-browser --ip=0.0.0.0 --LabApp.token="${TOKEN}"

CMD jupyter lab --port=8989 --no-browser --ip=0.0.0.0 --LabApp.token="${TOKEN}" --allow-root --NotebookApp.token=''
