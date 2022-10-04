#!/bin/bash -i
#https://www.ebi.ac.uk/ena/browser/view/PRJEB9691?show=reads
conda activate kingfisher

#already ran following line, redo if starting analysis from scratch
#kingfisher get -r -m ena-ascp #aws-http prefetch

#automating with run IDs once test download worked
while read line; do

	kingfisher get -r $line -m ena-ascp #aws-http prefetch

done < ids.txt
