#!/bin/bash -i
conda activate biom-env
#biom convert -i output-MG-combined-OTU-table/RAS_WGC_2016-2017.PacBio.merged.tsv -o output-MG-combined-OTU-table/RAS_WGC_2016-2017.PacBio.merged.biom --to-hdf5 --header-key taxonomy
biom convert -i output-MG-combined-OTU-table/RAS_WGC_2016-2017.PacBio.merged.SILVA-PR2-tax.tsv -o output-MG-combined-OTU-table/RAS_WGC_2016-2017.PacBio.merged.SILVA-PR2-tax.biom --to-hdf5 --header-key taxonomy

biom convert -i input-tag-ASV-tables/221129-1417_FRAM-Parada_2.0-fold-18S-correction_normalized_sequence_counts.tsv -o input-tag-ASV-tables/221129-1417_FRAM-Parada_2.0-fold-18S-correction_normalized_sequence_counts.biom --to-hdf5 --header-key taxonomy

biom convert -i input-tag-ASV-tables/221026-1554.FRAM-Wulf.18S.all-18S-seqs.with-SILVA-tax.tsv -o  input-tag-ASV-tables/221026-1554.FRAM-Wulf.18S.all-18S-seqs.with-SILVA-tax.biom --to-hdf5 --header-key taxonomy

