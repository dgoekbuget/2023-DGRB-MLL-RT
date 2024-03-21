#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate chipseq
#go to folder with bg files
for file in *.early.bg; do paste $file ${file%.early.bg}.late.bg | awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)}}' OFS='\t' > ${file%.early.bg}.RT.bg; done
echo -e "chr\tstart\tstop\t"`ls *RT.bg` | sed 's/\ /\t/g' > merge_RT.txt
bedtools unionbedg -filler "NA" -i *RT.bg >> merge_RT.txt

conda deactivate
