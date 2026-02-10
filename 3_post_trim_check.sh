#!/bin/bash -l

#SBATCH --time=2:00:00
#SBATCH --account=def-cdeboer
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --output=logs/3_post_trim_%j.out
#SBATCH --error=logs/3_post_trim_%j.err

module load fastqc

OUTDIR=outputs/fastqc/before_trim
DATA=outputs/fastqc/after_trim/trimmed_R1_paired.fq.gz
DATA2=outputs/fastqc/after_trim/trimmed_R2_paired.fq.gz
fastqc -o ${OUTDIR} -t 4 -f fastq ${DATA} ${DATA2}