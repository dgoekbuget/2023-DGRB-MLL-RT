#!/bin/bash
#SBATCH --mem=72G
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=12
#SBATCH --gres=scratch:200G

workdir="$1"
samplesheet="$workdir"/"$2"
wsize="$3"
windows=/blellochlab/data1/deniz/genomes/mm10."$wsize"kb.windows.bed
names=($(awk '{ print $1}' "$samplesheet" ))
name=${names[$SLURM_ARRAY_TASK_ID]}

source ~/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Map PE or SE data
x=`wc -l "$workdir"/"$name".bed | cut -d' ' -f 1`
bedtools intersect -sorted -c -b "$workdir"/"$name".bed -a "$windows" | awk -vx=$x '{print $1,$2,$3,$4*1e+06/x}' OFS='\t' > "$workdir"/"$name"."$wsize"kb.bg
rm "$workdir"/"$name".sam "$workdir"/"$name".bam

conda deactivate
