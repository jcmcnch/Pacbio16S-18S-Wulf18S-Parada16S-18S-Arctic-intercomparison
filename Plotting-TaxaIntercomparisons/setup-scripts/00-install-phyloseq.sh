#!/bin/bash -i

#set up a conda env as follows
conda create --name phyloseq-env
conda activate phyloseq-env
mamba install r-base r-tidyverse
