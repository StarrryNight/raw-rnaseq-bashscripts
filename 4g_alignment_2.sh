#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --account=rrg-cdeboer
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=40G
#SBATCH --output=logs/4c_alignment_%j.out
#SBATCH --error=logs/4c_alignment_%j.err

module load star/2.7.11b
module load samtools

# Define your paths
INDEX_DIR="data/hybrid_star_index"
FASTQ_R1="outputs/fastqc/after_trim/trimmed_R1_paired.fq.gz"
FASTQ_R2="outputs/fastqc/after_trim/trimmed_R2_paired.fq.gz"
OUT_PREFIX="outputs/alignment/yac_sample_1"

STAR --runThreadN 8 \
     --genomeDir $INDEX_DIR \
     --readFilesIn $FASTQ_R1 $FASTQ_R2 \
     --readFilesCommand zcat \
     --alignIntronMax 500000 \
     --outFileNamePrefix $OUT_PREFIX \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within  # This keeps unmapped reads in the BAM file for later