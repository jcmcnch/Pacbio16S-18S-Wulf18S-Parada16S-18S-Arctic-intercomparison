#!/bin/bash -i

mergedOutputWSC=WSC-merged-taxonomy
mergedOutputEGC=EGC-merged-taxonomy
mkdir -p $mergedOutputWSC
mkdir -p $mergedOutputEGC

for inputFolder in EGC-classified-SILVA138.1 WSC-classified-SILVA138.1 ; do

        inputSILVA=$inputFolder
	inputPR2=`echo $inputSILVA | sed 's/SILVA138.1/PR2-4.14.0/'`
	stationName=`echo $inputFolder | cut -f1 -d-`

		#loop over SILVA taxonomy, getting IDs that need to be extracted from PR2 to replace SILVA taxonomy. Same overall steps as in ASV tag pipeline to make it methodologically consistent.
                for item in `ls $inputSILVA/*.classified.SILVA138.1.tsv`; do

			SampleID=`basename $item _SSU_rRNA.classified.SILVA138.1.tsv`
			
			#collect IDs for subsequent grepping
			IDs4PR2=`basename $item .tsv`.IDs4PR2.tsv
			grep "d:Eukaryota" $item | cut -f1 > $IDs4PR2
			grep "o:Chloroplast" $item | cut -f1 >> $IDs4PR2

			#specify output file
			mergedOutput=$stationName-merged-taxonomy/$SampleID.merged-tax.tsv

			#take taxonomy from SILVA first using grep -v
			grep -v -f $IDs4PR2 $item > $mergedOutput

			#now get Chloroplast / Euk taxonomy from PR2
			grep -f $IDs4PR2 $inputPR2/${SampleID}_SSU_rRNA.classified.PR2-4.14.0.tsv >> $mergedOutput

			#remove unneeded IDs file
			rm -f $IDs4PR2

        done
done

