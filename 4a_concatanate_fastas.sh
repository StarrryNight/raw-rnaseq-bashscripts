#!/bin/bash -l

#SBATCH --time=2:00:00
#SBATCH --account=def-cdeboer
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --output=4b_concatanate_fastas_%j.out
#SBATCH --error=4b_concatanate_fastas_%j.err

HUMAN_FA=data/hg38.fa
YEAST_FA=data/sacCer3.fa

(cat ${YEAST_FA}; echo ""; cat ${HUMAN_FA}) > data/yeast_human_hybrid.fa