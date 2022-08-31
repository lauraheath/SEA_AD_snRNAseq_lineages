FROM rocker/rstudio:4.1.0

RUN apt-get update -y\
&& apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev zlib1g-dev libbz2-dev liblzma-dev build-essential libglpk40\
&& apt-get install -y curl libcurl4-openssl-dev libxml2-dev libcairo2-dev libxt-dev r-cran-cairo libgsl-dev\
&& apt-get install -y git openssh-client nano less libxt6 libcurl4-openssl-dev python3 python3-pip

RUN pip3 install synapseclient

RUN R -e "install.packages('remotes')" \
&& R -e "install.packages('BiocManager')"\
&& R -e "install.packages('devtools')"\
&& R -e "install.packages('R.utils')"\
&& R -e "BiocManager::install('biomaRt')"\
&& R -e "BiocManager::install('monocle')"\
&& R -e "install.packages('tidyverse')"\
&& R -e "install.packages('caret')"\
&& R -e "install.packages('data.table')"\
&& R -e "install.packages('ggplot2')"\
&& R -e "install.packages('synapser', repos = c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"\
&& R -e "devtools::install_github('lhe17/nebula')"\
&& R -e "devtools::install_github('satijalab/seurat-wrappers')"\
&& R -e "install.packages('knitr')"\
&& R -e "BiocManager::install('limma')"\
&& R -e "install.packages('Matrix')"\
&& R -e "install.packages('patchwork')"\
&& R -e "install.packages('RColorBrewer')"\
&& R -e "install.packages('Seurat')"\
&& R -e "devtools::install_github('cole-trapnell-lab/leidenbase')"
