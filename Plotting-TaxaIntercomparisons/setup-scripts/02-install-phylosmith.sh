#!/bin/bash -i
#instructions for installing phylosmith
#tried manual "R" way, but failed. So resorted to mamba/conda
sudo apt install libmysqlclient-dev libgdal-dev libudunits2-dev
mamba install -c conda-forge r-devtools
mamba install -c conda-forge r-rcppeigen
mamba install -c conda-forge r-rcppparallel
mamba install -c conda-forge r-rtsne
mamba install -c conda-forge r-ggforce
mamba install -c conda-forge r-units
mamba install -c conda-forge r-sf

#then run following:
remotes::install_github('schuyler-smith/phylosmith') 
