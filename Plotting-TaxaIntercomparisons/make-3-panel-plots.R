library(phyloseq)
library(tidyverse)
library(phylosmith)
library(gridExtra)
#import files from biom, which works now
MG <- import_biom("output-MG-combined-OTU-table/RAS_WGC_2016-2017.PacBio.merged.biom")
colnames(tax_table(MG)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
Parada <- import_biom("input-tag-ASV-tables/221129-1417_FRAM-Parada_2.0-fold-18S-correction_normalized_sequence_counts.biom")
#note the taxonomy is different here because of PR2 having 8 levels
colnames(tax_table(Parada)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Subspecies")
Wulf <- import_biom("input-tag-ASV-tables/230307-2053.FRAM-Wulf.18S.all-18S-seqs.with-SILVA-tax.biom")
colnames(tax_table(Wulf)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#Wulf <- subset_taxa(Wulf, Domain != "unassigned")

#next import sample data
MGmeta <- read.delim("metadata/221110_sample_metadata.PacBio.tsv", sep = "\t", header = TRUE)
rownames(MGmeta) <- MGmeta$SampleID

ParadaMeta <- read.delim("metadata/221110_sample_metadata.Parada.tsv", sep = "\t", header = TRUE)
rownames(ParadaMeta) <- ParadaMeta$SampleID

WulfMeta <- read.delim("metadata/221110_sample_metadata.Wulf.tsv", sep = "\t", header = TRUE)
rownames(WulfMeta) <- WulfMeta$SampleID

#and delete the first column because it is now redundant
MGmeta <- MGmeta %>% 
  select(-SampleID)
ParadaMeta <- ParadaMeta %>% 
  select(-SampleID)
WulfMeta <- WulfMeta %>% 
  select(-SampleID)

#import into phyloseq object, then merge
MGmeta <- sample_data(MGmeta)
ParadaMeta <- sample_data(ParadaMeta)
MG <- merge_phyloseq(MG, MGmeta)
Parada <-merge_phyloseq(Parada, ParadaMeta)
WulfMeta <- sample_data(WulfMeta)
Wulf <- merge_phyloseq(Wulf, WulfMeta)

#manually set sample order by date
MG <- set_sample_order(MG, c('RAS_WSC_09_2016_1', 'RAS_WSC_10_2016', 'RAS_WSC_11_2016', 'RAS_WSC_12_2016_1', 'RAS_WSC_03_2017_2', 'RAS_WSC_05_2017_1', 'RAS_WSC_06_2017_1', 'RAS_WSC_07_2017_1', 'RAS_WSC_07_2017_2', 'RAS_WSC_07_2017_3'))

#no need to set by date since ordered in metadata. Couldn't get order function to work for some unknown reason. Probably to do with the smaller number of samples.
Parada <- subset_samples(Parada, Sample_type=="F4")
Wulf <- subset_samples(Wulf, Sample_type=="F4")
Wulf <- prune_taxa(taxa_sums(Wulf) > 0, Wulf)

#then do some plotting
MGline <- abundance_lines(MG, classification = 'Domain', relative_abundance = TRUE, treatment = "Sample_type", sample_labels = sample_data(MG)$date) + theme(axis.text.x = element_blank()) + scale_y_continuous(trans = "log10")
ParadaLine <- abundance_lines(Parada, classification = 'Domain', relative_abundance = TRUE, treatment = "Sample_type", sample_labels = sample_data(Parada)$date) + theme(axis.text.x = element_blank()) + scale_y_continuous(trans = "log10")
WulfLine <- abundance_lines(Wulf, classification = 'Domain', relative_abundance = TRUE, treatment = "Sample_type", sample_labels = sample_data(Wulf)$date) #+ theme(axis.text.x = element_blank())
grid.arrange(MGline, ParadaLine, ncol=1)

#TODO:
#Enforce same width
#Unify taxonomy so that colours can be unified and tax_glom etc will be applied similarly
#Play with different ways of plotting axes so that lower abundance taxa are visible
#Think about transformation of data as suggested by Matthias

#Still can't figure out why the otu_table slot is empty. This *might* be a phylosmith error. I can try to open an issue there, because googling about isn't really helping me.

phylogeny_profile(Wulf, classification = 'Domain', treatment = "Sample_type", subset = NULL, merge = TRUE, relative_abundance = TRUE, colors = 'default')
WulfSubset <- subset_taxa(Wulf, Domain=="Eukaryota")
plot_bar(Wulf.Dinoflagellata,"date","Genus")
#MGmeta <- MGmeta[order(as.Date(MGmeta$date, format="%m/%d/%Y")),]
#ParadaMeta <- ParadaMeta[order(as.Date(ParadaMeta$date, format="%m/%d/%Y")),]
#Parada <- set_sample_order(Parada, c('RAS_WSC_09_2016_1', 'RAS_WSC_10_2016', 'RAS_WSC_11_2016', 'RAS_WSC_12_2016_1', 'RAS_WSC_03_2017_2', 'RAS_WSC_05_2017_1', 'RAS_WSC_06_2017_1', 'RAS_WSC_07_2017_1', 'RAS_WSC_07_2017_2', 'RAS_WSC_07_2017_3'))
#Parada <- set_sample_order(Parada, c("3409-FRAM-bactV4V5-F4-S-1-5plus6-S75-L001","3411-FRAM-bactV4V5-F4-S-1-9plus10-S77-L001","3412-FRAM-bactV4V5-F4-S-1-11plus12-S78-L001","3413-FRAM-bactV4V5-F4-S-1-13plus14-S79-L001","3418-FRAM-bactV4V5-F4-S-1-23plus24-S84-L001","3423-FRAM-bactV4V5-F4-S-1-33plus34-S89-L001","3425-FRAM-bactV4V5-F4-S-1-37plus38-S91-L001","3427-FRAM-bactV4V5-F4-S-1-41plus42-S93-L001","3428-FRAM-bactV4V5-F4-S-1-43plus44-S94-L001","3429-FRAM-bactV4V5-F4-S-1-45plus46-S95-L001"))


#export to tsv to troubleshoot
write_biom_tsv <- function(ps, file, sep = "; ") {
  phyloseq::otu_table(ps) %>%
    as.data.frame() %>%
    rownames_to_column("#OTU ID") %>%
    left_join(phyloseq::tax_table(ps) %>%
                as.data.frame() %>%
                rownames_to_column("#OTU ID") %>%
                tidyr::unite("taxonomy", !`#OTU ID`, sep = sep)) -> phyloseq_biom
  
  write_tsv(phyloseq_biom, file = file)
}

write_biom_tsv(Wulf, "Wulf.test.tsv", sep = "; ")
