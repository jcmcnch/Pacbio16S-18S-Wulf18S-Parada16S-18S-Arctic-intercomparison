#Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
sed -re 's/\([0-9]{1}.[0-9]{2}\)//g' $1 | cut -f2 | sort | uniq -c  | awk '{print $2,"\t",$1}' | sed 's/^ \t/NA\t/g'

#                "tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.matches} ; " #Process output into tsv format to stdout
#sed -re 's/\([0-9]{{1}}\.[0-9]{{2}}\)//g' {input.mismatches} | tee {output.taxTableMismatches} |" #Remove confidence estimations from VSEARCH output, keep a copy for later steps but also pipe to subsequent commands
#                "cut -f2 | sort | cut -d, -f1-4 | sort | uniq -c | " #Take only tax column, collapse to order level, then count unique occurrences
#                "tail -f -n +2 | awk '{{print $1,\"\t\",$2}}' > {output.mismatches}
