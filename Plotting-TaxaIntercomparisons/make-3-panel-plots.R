library(phyloseq)
library(tidyverse)
library(phylosmith)
#import files from biom, which works now
MG <- import_biom("output-MG-combined-OTU-table/RAS_WGC_2016-2017.PacBio.merged.biom")
colnames(tax_table(MG)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
Parada <- import_biom("input-tag-ASV-tables/221129-1417_FRAM-Parada_2.0-fold-18S-correction_normalized_sequence_counts.biom")
Wulf <- import_biom("input-tag-ASV-tables/221026-1554.FRAM-Wulf.18S.all-18S-seqs.with-SILVA-tax.biom")

#TODO:
#change taxonomy headers from "Rank1, Rank2, ... etc" to something more biologically realistic.

#next import sample data
MGmeta <- read.delim("metadata/221110_sample_metadata.PacBio.tsv", sep = "\t", header = TRUE)
rownames(MGmeta) <- MGmeta$SampleID

#and delete the first column because it is now redundant
MGmeta <- MGmeta %>% 
  select(-SampleID)

MGmeta <- MGmeta[order(as.Date(MGmeta$date, format="%m/%d/%Y")),]

#import into phyloseq object, then merge
MGmeta <- sample_data(MGmeta)
MG <- merge_phyloseq(MG, MGmeta)
MG <- set_sample_order(MG, c('RAS_WSC_09_2016_1', 'RAS_WSC_10_2016', 'RAS_WSC_11_2016', 'RAS_WSC_12_2016_1', 'RAS_WSC_03_2017_2', 'RAS_WSC_05_2017_1', 'RAS_WSC_06_2017_1', 'RAS_WSC_07_2017_1', 'RAS_WSC_07_2017_2', 'RAS_WSC_07_2017_3'))

#then do some plotting
abundance_lines(MG, classification = 'Domain', relative_abundance = TRUE, treatment = "Sample_type", sample_labels = sample_data(MG)$date)
phylogeny_profile(MG, classification = 'Domain', treatment = "Sample_type", subset = NULL, merge = TRUE, relative_abundance = TRUE, colors = 'default')
