#!/bin/bash -i

for inputFolder in EGC-merged-taxonomy WSC-merged-taxonomy ; do

        outdir=$inputFolder

                for item in `ls $inputFolder/*.merged-tax.tsv`; do

                        SampleID=`basename $item .merged-tax.tsv`
                        outfile=${SampleID}_SSU_rRNA.merged-tax.taxtable.tsv
                        printf "SampleID\t$SampleID\n" > $inputFolder/$outfile
                        sed -re 's/\([0-9]{1}.[0-9]{2}\)//g' $item | cut -f2 | cut -f1-5 -d, | sort | uniq -c | awk '{print $2,"\t",$1}' | sed 's/^ \t/d:Unclassified\t/g' >> $inputFolder/$outfile

        done
done

#for inputFolder in EGC-classified-SILVA138.1 WSC-classified-SILVA138.1 ; do

#        outdir=$inputFolder

#                for item in `ls $inputFolder/*.classified.SILVA138.1.tsv`; do

#			SampleID=`basename $item _SSU_rRNA.classified.SILVA138.1.tsv`
#			outfile=`basename $item .tsv`.taxtable.tsv
#			printf "SampleID\t$SampleID\n" > $inputFolder/$outfile
#			sed -re 's/\([0-9]{1}.[0-9]{2}\)//g' $item | cut -f2 | cut -f1-5 -d, | sort | uniq -c | awk '{print $2,"\t",$1}' | sed 's/^ \t/d:Unclassified\t/g' >> $inputFolder/$outfile

#        done
#done


#Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
#sed -re 's/\([0-9]{1}.[0-9]{2}\)//g' $1 | cut -f2 | sort | uniq -c  | awk '{print $2,"\t",$1}' | sed 's/^ \t/NA\t/g'

#                "tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.matches} ; " #Process output into tsv format to stdout
#sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.mismatches} | tee {output.taxTableMismatches} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
#                "cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
#                "tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.mismatches}
