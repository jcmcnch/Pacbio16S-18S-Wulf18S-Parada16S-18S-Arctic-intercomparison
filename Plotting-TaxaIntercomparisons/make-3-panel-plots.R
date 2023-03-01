library(phyloseq)
library(tidyverse)
#import file
MG <- import_biom("output-MG-combined-OTU-table/RAS_WGC_2016-2017.PacBio.merged.biom")
Parada <- import_biom("input-tag-ASV-tables/221129-1417_FRAM-Parada_2.0-fold-18S-correction_normalized_sequence_counts.biom")
Wulf <- import_biom("input-tag-ASV-tables/221026-1554.FRAM-Wulf.18S.all-18S-seqs.with-SILVA-tax.biom")

