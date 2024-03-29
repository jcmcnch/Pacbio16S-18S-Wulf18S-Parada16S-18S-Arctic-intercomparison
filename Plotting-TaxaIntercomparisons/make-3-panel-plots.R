library(phyloseq)
library(tidyverse)
library(phylosmith)
library(gridExtra)
library(gtable)
library(grid)
#import files from biom, which works now
MG <- import_biom("output-MG-combined-OTU-table/RAS_WGC_2016-2017.PacBio.merged.SILVA-PR2-tax.biom")
#colnames(tax_table(MG)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
Parada <- import_biom("input-tag-ASV-tables/221129-1417_FRAM-Parada_2.0-fold-18S-correction_normalized_sequence_counts.biom")
#note the taxonomy is different here because of PR2 having 8 levels
#colnames(tax_table(Parada)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species","Subspecies")
Wulf <- import_biom("input-tag-ASV-tables/230307-2053.FRAM-Wulf.18S.all-18S-seqs.with-PR2-tax.biom")
#colnames(tax_table(Wulf)) <- c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species", "Subspecies")

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
Wulf <- set_sample_order(Wulf, c('ERR5533347','ERR5533345','ERR5533344','ERR5533343','ERR5533338','ERR5533333','ERR5533331','ERR5533329','ERR5533328','ERR5533327'))

Wulf = relative_abundance(Wulf)
Parada = relative_abundance(Parada)
MG = relative_abundance(MG)

taxa_of_interest = c("c:Bacillariophyta","c:Bacillariophyta:plas","c:Arthropoda") #"c:Spirotrichea",

MG <- conglomerate_taxa(MG, "Rank4")
MG_subset <- taxa_extract(MG, taxa_of_interest)
MG_subset <- microbiome::transform(MG_subset,"clr")

taxa_of_interest_Parada = c("Bacillariophyta","Bacillariophyta:plas","Arthropoda") #"Spirotrichea",
Parada_subset <- taxa_extract(Parada, taxa_of_interest_Parada)
Parada_subset <- microbiome::transform(Parada_subset,"clr")

taxa_of_interest_Wulf = c("Bacillariophyta","Bacillariophyta:plas","Arthropoda") #"Spirotrichea",
Wulf_subset <- taxa_extract(Wulf, taxa_of_interest_Wulf)
Wulf_subset <- microbiome::transform(Wulf_subset, "clr")

#MG_subset <- subset_taxa(MG, Rank4=="c:Bacillariophyta" | Rank3=="p:Ciliophora" | Rank5=="o:Crustacea")
#Parada_subset <- subset_taxa(Parada, Rank5=="Bacillariophyta_X:plas" | Rank3=="Ciliophora" | Rank5=="Crustacea")
#Wulf_subset <- subset_taxa(Wulf, Rank4=="Bacillariophyta" | Rank3=="Ciliophora"  | Rank5=="Crustacea")

#PacBio
#make plot and set theme
MGline <- abundance_lines(MG_subset, classification = 'Rank4', relative_abundance = FALSE, treatment = "method", 
                          sample_labels = sample_data(MG)$date) + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), legend.position = "none") 
  #+ scale_y_continuous(trans = "log10")

#apply custom colors
color_map <- data.frame(taxa = taxa_of_interest,
                        colors = phylosmith:::create_palette(length(taxa_of_interest)))
MG_colors <- color_map$colors[match(unique(MGline$data$OTU), color_map$taxa, FALSE)]
MGline <- MGline + ggplot2::scale_color_manual(values = MG_colors) 

#Parada primers
ParadaLine <- abundance_lines(Parada_subset, classification = 'Rank4', relative_abundance = FALSE, treatment = "method", 
                              sample_labels = sample_data(Parada)$date) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.title = element_blank()) 
  #+ scale_y_continuous(trans = "log10")
color_map_Parada <- data.frame(taxa = taxa_of_interest_Parada,
                        colors = phylosmith:::create_palette(length(taxa_of_interest_Parada)))
Parada_colors <- color_map_Parada$colors[match(unique(ParadaLine$data$Rank4), color_map_Parada$taxa, FALSE)]
ParadaLine <- ParadaLine + ggplot2:::scale_color_manual(values = Parada_colors)

#Wulf 18S primers
#cannot get rid of legend here for some odd reason. Works with MG, but not wulf. Maybe a weird typo somewhere.
WulfLine <- abundance_lines(Wulf_subset, classification = 'Rank4', relative_abundance = FALSE, treatment = "method", 
                            sample_labels = sample_data(Wulf)$date) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.title = element_blank())#, legend.position = "none") 
#+ scale_y_continuous(trans = "log10")
color_map_Wulf <- data.frame(taxa = taxa_of_interest_Wulf,
                        colors = phylosmith:::create_palette(length(taxa_of_interest_Wulf)))
Wulf_colors <- color_map_Wulf$colors[match(unique(WulfLine$data$Rank4), color_map_Wulf$taxa, FALSE)]
WulfLine <- WulfLine + ggplot2:::scale_color_manual(values = Wulf_colors)
g1 <- ggplotGrob(MGline)
g2 <- ggplotGrob(ParadaLine)
g3 <- ggplotGrob(WulfLine)
g1$widths=g2$widths=g3$widths
g <- rbind(g1, g2, g3, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
g$heights <- unit.pmax(g1$heights,g2$heights,g3$heights)
grid.newpage()
grid.arrange(g1, g2, g3, ncol=1)