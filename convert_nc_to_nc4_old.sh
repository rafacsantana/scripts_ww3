#!/bin/bash

#install
#brew install ghostscript


FILES=`ls roms*.nc`
for f in ${FILES[@]}; do

 echo "Processing file: ${f}"

 base=`basename ${f} .nc`
 #base=small_$base

 #convert ${f} ${base}.eps
 #convert ${f} -resize 50% ${base}.eps

 #compresspdf "IMG_7129.pdf" "IMG_7129r.pdf" printer
 #compresspdf ${f} ${base}.pdf printer

 #nccopy -d1 classic.nc compressed.nc 
 nccopy -d1 ${f} ${base}.nc4 
 mv  ${base}.nc4 ${f} 

done

