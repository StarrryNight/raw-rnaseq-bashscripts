#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --account=rrg-cdeboer
#SBATCH --job-name=bamcoverage
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --output=logs/6_match_fa_%j.out
#SBATCH --error=logs/6_match_fa_%j.err

REF_FASTA="data/hg38.fa"

read TOP_CHR START END < my_yac.bed

module load samtools
# Extract the sequence
samtools faidx ${REF_FASTA} ${TOP_CHR}:${START}-${END} > final_output/human_yac_insert.fa