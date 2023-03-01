#!/bin/bash -i
conda activate phyloseq-env

cd input-MG-sample-tables/WSC/ 
~/github/eASV-pipeline-for-515Y-926R/qiime2-2022.2-DADA2-SILVA138.1-PR2_4.14.0/02-utility-scripts/merge-SINTAX-tsv-phyloflash-and-output-as-OTU-table.R --inputTSV "RAS_WSC_03_2017_2_SSU_rRNA.classified.SILVA138.1.taxtable.tsv,RAS_WSC_05_2017_1_SSU_rRNA.classified.SILVA138.1.taxtable.tsv,RAS_WSC_06_2017_1_SSU_rRNA.classified.SILVA138.1.taxtable.tsv,RAS_WSC_07_2017_1_SSU_rRNA.classified.SILVA138.1.taxtable.tsv,RAS_WSC_07_2017_2_SSU_rRNA.classified.SILVA138.1.taxtable.tsv,RAS_WSC_07_2017_3_SSU_rRNA.classified.SILVA138.1.taxtable.tsv,RAS_WSC_09_2016_1_SSU_rRNA.classified.SILVA138.1.taxtable.tsv,RAS_WSC_10_2016_SSU_rRNA.classified.SILVA138.1.taxtable.tsv,RAS_WSC_11_2016_SSU_rRNA.classified.SILVA138.1.taxtable.tsv,RAS_WSC_12_2016_1_SSU_rRNA.classified.SILVA138.1.taxtable.tsv" --OTUout ../../output-MG-combined-OTU-table/RAS_WGC_2016-2017.PacBio.merged.tsv
