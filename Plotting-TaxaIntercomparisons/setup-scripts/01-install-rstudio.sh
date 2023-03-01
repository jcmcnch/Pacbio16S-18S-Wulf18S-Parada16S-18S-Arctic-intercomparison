#!/bin/bash -i

conda activate phyloseq-env
#install rstudio
mamba install -c conda-forge rstudio-desktop
mamba install -c r r-codetools
mamba install -c conda-forge r-rcpp
mamba install -c conda-forge r-rcurl
