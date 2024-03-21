#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ucsc
for i in *Loess.bg; do echo $i && bedGraphToBigWig $i /blellochlab/data1/deniz/genomes/mm10.sorted.chrom.sizes ${i%.bg}.bw; done
conda deactivate
mkdir -p bw
mv *.bw bw/
