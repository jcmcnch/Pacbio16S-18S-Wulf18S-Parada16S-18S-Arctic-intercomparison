#!/bin/bash -i

conda activate qiime2-2019.4

mkdir 01-trimmed
mkdir -p logs/01-trimmed/

for item in `ls 00-raw/*_1.fastq.gz` 

	do

	filestem=`basename $item _1.fastq.gz`
	R1=00-raw/${filestem}_1.fastq.gz
	R2=00-raw/${filestem}_2.fastq.gz

        cutadapt --no-indels --pair-filter=any --error-rate=0.2 --discard-untrimmed \
	-g ^GCGGTAATTCCAGCTCCAA -G ^ACTTTCGTTCTTGATYRR \
	-o 01-trimmed/${filestem}.R1.trimmed.fastq \
	-p 01-trimmed/${filestem}.R2.trimmed.fastq $R1 $R2 \
	2>&1 | tee -a logs/01-trimmed/${filestem}.cutadapt.stderrout

done

conda deactivate
