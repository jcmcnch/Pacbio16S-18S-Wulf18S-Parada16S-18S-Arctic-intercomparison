library(phyloseq)
library(tidyverse)
#import file
MG <- import_biom("output-MG-combined-OTU-table/RAS_WGC_2016-2017.PacBio.merged.biom")
Parada <- import_biom("foo.biom")