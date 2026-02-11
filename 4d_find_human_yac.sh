#!/bin/bash -l

#SBATCH --time=2:00:00
#SBATCH --account=rrg-cdeboer
#SBATCH --job-name=find-human_yac
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=40G
#SBATCH --output=logs/4d_find_human_yac_%j.out
#SBATCH --error=logs/4d_find_human_yac_%j.err

module load samtools

YAC_SAMPLE_FILE=outputs/alignment/yac_sample_1Aligned.sortedByCoord.out.bam

TOP_CHR=$(samtools idxstats ${YAC_SAMPLE_FILE} | \
grep "^chr" | \
awk '$2 > 10000000' | \
sort -k3,3rn | \
head -n 1 |\
cut -f1)

echo "Top chromosome is: $TOP_CHR"
#TODO rename duplicated chrom names
# This sorts all 1.2M positions and looks at the middle 
# coordinates to find the 'core' of the YAC.
samtools view ${YAC_SAMPLE_FILE} "$TOP_CHR" | \
    awk '{print $4}' | \
    sort -n | \
    sed -n '100000,1100000p' | \
    sed -n '1p;$p'