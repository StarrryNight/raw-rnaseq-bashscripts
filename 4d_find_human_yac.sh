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

samtools index ${YAC_SAMPLE_FILE}
samtools idxstats ${YAC_SAMPLE_FILE} | \
grep "^chr" | \
awk '$2 > 10000000' | \
sort -k3,3rn | \
head -n 5
