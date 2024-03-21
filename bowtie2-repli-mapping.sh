#!/bin/bash
#SBATCH --mem=72G
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=12
#SBATCH --gres=scratch:200G

workdir="$1"
samplesheet="$workdir"/"$2"
bt2index="/blellochlab/data1/deniz/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
windows="/blellochlab/data1/deniz/genomes/mm10.50kb.windows.bed"
names=($(awk '{ print $1}' "$samplesheet" ))
read1s=($(awk '{ print $2}' "$samplesheet" ))
read2s=($(awk '{ print $3}' "$samplesheet" ))
name=${names[$SLURM_ARRAY_TASK_ID]}
read1=${read1s[$SLURM_ARRAY_TASK_ID]}
read2=${read2s[$SLURM_ARRAY_TASK_ID]}

source ~/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Map PE or SE data
if [ "$read2" == "SE" ]; then
  echo "$name was sequenced 'SE'"
  bowtie2 -x "$bt2index" -p 12 --no-mixed --no-discordant --reorder -U $read1 -S "$workdir"/"$name".sam 2>> "$workdir"/"$name"_mapping_log.txt
else
  echo "$name was sequenced PE"
  bowtie2 -x "$bt2index" -p 12 --no-mixed --no-discordant --reorder -X 1000 -1 $read1 -2 $read2 -S "$workdir"/"$name".sam 2>> "$workdir"/"$name"_mapping_log.txt
fi

samtools view -bSq 20 "$workdir"/"$name".sam > "$workdir"/"$name".bam
samtools sort -o "$workdir"/"$name"_srt.bam "$workdir"/"$name".bam
samtools rmdup -S "$workdir"/"$name"_srt.bam "$workdir"/"$name"_rmdup.bam
bamToBed -i "$workdir"/"$name"_rmdup.bam | cut -f 1,2,3,4,5,6 | sort -T . -k1,1 -k2,2n -S 5G > "$workdir"/"$name".bed
x=`wc -l "$workdir"/"$name".bed | cut -d' ' -f 1`
bedtools intersect -sorted -c -b "$workdir"/"$name".bed -a "$windows" | awk -vx=$x '{print $1,$2,$3,$4*1e+06/x}' OFS='\t' > "$workdir"/"$name".bg
rm "$workdir"/"$name".sam "$workdir"/"$name".bam

conda deactivate
