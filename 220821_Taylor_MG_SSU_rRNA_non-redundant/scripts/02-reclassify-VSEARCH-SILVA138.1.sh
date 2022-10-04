#!/bin/bash -i
conda activate vsearch-env

for inputFolder in EGC WSC ; do 

	mkdir -p $inputFolder-classified-SILVA138.1
	outdir=$inputFolder-classified-SILVA138.1

		for item in `ls $inputFolder/*fa`; do

			outfile=`basename $item .fa`.classified.SILVA138.1.tsv
			vsearch --sintax $item \
			       	--db /home/db/VSEARCH/220830_SILVA138.1_VSEARCH-formatted.udb \
			        --tabbedout $outdir/$outfile --threads 1 --sintax_cutoff 0 \
				--top_hits_only --topn 1 --notrunclabels

	done
done
