#!/bin/bash -i

#set up a conda env as follows
conda create --name phyloseq-env
conda activate phyloseq-env
mamba install r-base r-tidyverse

#entered interactive R session
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.15")
#BiocManager::install("phyloseq")
