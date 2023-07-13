library(phyloseq)
library(tidyverse)
library(phylosmith)
library(gridExtra)
#import files from biom, which works now
Wulf <- import_biom("input-tag-ASV-tables/230307-2053.FRAM-Wulf.18S.all-18S-seqs.with-SILVA-tax.biom")
colnames(tax_table(Wulf)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
Wulf <- subset_taxa(Wulf, Domain != "unassigned")

#next import sample data
WulfMeta <- read.delim("metadata/221110_sample_metadata.Wulf.tsv", sep = "\t", header = TRUE)
rownames(WulfMeta) <- WulfMeta$SampleID

#and delete the first column because it is now redundant
WulfMeta <- WulfMeta %>% 
  select(-SampleID)

#import into phyloseq object, then merge
WulfMeta <- sample_data(WulfMeta)
Wulf <- merge_phyloseq(Wulf, WulfMeta)

#no need to set by date since ordered in metadata. Couldn't get order function to work for some unknown reason. Probably to do with the smaller number of samples.
Wulf <- subset_samples(Wulf, Sample_type=="F4")
Wulf <- prune_taxa(taxa_sums(Wulf) > 0, Wulf)
Wulf <- set_sample_order(Wulf, c('ERR5533347','ERR5533345','ERR5533344','ERR5533343','ERR5533338','ERR5533333','ERR5533331','ERR5533329','ERR5533328','ERR5533327'))
Wulf <- relative_abundance(Wulf)

#then do some plotting
plot_bar(Wulf, "date", "Abundance", "Domain")
taxa_core_graph(Wulf)
taxa_abundance_bars(Wulf, classification = "Phylum")

abundance_lines(Wulf, treatment = "Sample_type", sample_labels = sample_data(Wulf)$date) + theme(axis.text.x = element_blank())
#abundance_heatmap(Wulf)
