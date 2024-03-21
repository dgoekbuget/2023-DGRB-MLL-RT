#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate deeptools
#go to folder with bw files

for i in *.n1.*Loess.bw; do
	bigwigCompare -b1 $i -b2 ${i%.n1.RT.Loess.bw}.n2.RT.Loess.bw --operation add -p 12 -o $i.temp -of bigwig -bs 5000
	bigwigCompare -b1 $i.temp -b2 ${i%.n1.RT.Loess.bw}.n3.RT.Loess.bw --operation mean -p 12 -o ${i%.n1.RT.Loess.bw}.mean.RT.Loess.bw -of bigwig -bs 5000 --scaleFactors 0.5:1
done

rm *.temp

conda deactivate
