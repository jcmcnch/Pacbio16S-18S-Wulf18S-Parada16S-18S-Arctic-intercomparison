#!/bin/bash -i
conda activate vsearch-env

for inputFolder in EGC WSC ; do 

	mkdir -p $inputFolder-classified-PR2-4.14.0
	outdir=$inputFolder-classified-PR2-4.14.0

		for item in `ls $inputFolder/*fa`; do

			outfile=`basename $item .fa`.classified.PR2-4.14.0.tsv
			vsearch --sintax $item \
			       	--db /home/db/VSEARCH/220830_PR2-4.14.0_VSEARCH-formatted.udb \
			        --tabbedout $outdir/$outfile --threads 1 --sintax_cutoff 0 \
				--top_hits_only --topn 1 --notrunclabels

	done
done
