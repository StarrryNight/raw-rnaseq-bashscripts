#!/bin/bash -l

#SBATCH --time=2:00:00
#SBATCH --account=def-cdeboer
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --output=logs/fastqc_raw_%j.out
#SBATCH --error=logs/fastqc_raw_%j.err

module load fastqc

OUTDIR=outputs/fastqc/before_trim
DATA=data/1951f6-DL-8_10_S13_L007_R1_001.fastq.gz

fastqc -o ${OUTDIR} -t 4 -f fastq ${DATA}