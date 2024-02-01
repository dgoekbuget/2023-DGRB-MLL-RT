#!/bin/bash
#SBATCH --mem=72G
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=12

workdir="$1"
files=($(awk '{print $1}' "$workdir"/"$2"))
reference="$workdir"/"$3"
genome="/blellochlab/data1/deniz/genomes/mm10.20kb.windows.bed"
outdir="$workdir"/background

mkdir -p "$outdir"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

mark="${files[$SLURM_ARRAY_TASK_ID]}"
echo "handling $mark"

if find /blellochlab/data1/Boileau/data/cutandtag/naivetoform/*/03_peak_calling/05_consensus_peaks -name *aive*"$mark"*.consensus.peaks.bed | while read filename; do cat "$filename" | bedtools subtract -a "$genome" -b stdin | awk -v OFS="\t" '{print NR,$1,$2,$3,"-"}' > "$outdir"/"$mark".temp.saf; done; then
echo "Background reference file for $mark generated"
fi 

{ echo -e "GeneID\tChr\tStart\tEnd\tStrand"; cat "$outdir"/"$mark".temp.saf; } > "$outdir"/"$mark".background.saf

rm "$outdir"/"$mark".temp.saf

conda deactivate
conda activate counts

bamfiles=($(find /blellochlab/data1/Boileau/data/cutandtag/naivetoform/all_bam_files -name *aive*"$mark"*.target.markdup.bam))
echo "bamfiles are: ${bamfiles[@]}"

if featureCounts -F 'SAF' -s 0 -M -p --fraction -a "$outdir"/"$mark".background.saf -o "$outdir"/"$mark".background.counts.txt ${bamfiles[@]}; then echo "Count file for $mark generated!"
fi

blength=$(awk '{NR > 2;sum+=$6} END{print sum;}' "$outdir"/"$mark".background.counts.txt)
sum=$(awk '{NR > 2;sum+=$7} END{print sum;}' "$outdir"/"$mark".background.counts.txt)
sum2=$(awk '{NR > 2;sum+=$8} END{print sum;}' "$outdir"/"$mark".background.counts.txt)
background=$(echo "scale=10 ; ($sum/2+$sum2/2)/$blength" | bc)
plength=($(find /blellochlab/data1/Boileau/data/cutandtag/naivetoform/*/03_peak_calling/05_consensus_peaks -name *aive*"$mark"*.consensus.peaks.bed | while read filename; do cat "$filename" | awk -v OFS='\t' '$10 != 1 {print}' | awk '{print int($3-$2)}'; done))

echo "Background length is: $blength, Sum R1 is: $sum, Sum R2 is: $sum2, Background for $mark is: $background"
echo "Peak lengths for $mark starts with $plength"

conda deactivate
conda activate chipseq

if find /blellochlab/data1/Boileau/data/cutandtag/naivetoform/*/03_peak_calling/05_consensus_peaks -name *aive*"$mark"*.consensus.peaks.bed | while read filename; do cat "$filename" | awk -v OFS='\t' '$10 != 1 {print}'| awk -v OFS='\t' -F'[\t,]' '{print $1,$2,$3,$8,$9}' | paste - <(echo ${plength[@]}) | awk -v i="$background" -v OFS='\t' '{print $1,$2,$3,$4-i*$6,$5-i*$6}' | bedtools map -a $reference -b stdin -c 4,5 -o sum -F 0.5 -null 0 > "$workdir"/"$mark".RTdomain.backgroundCorSum.bed; done
then echo "Background corrected peak sum bed file for $mark generated"
fi

conda deactivate
#EOF
