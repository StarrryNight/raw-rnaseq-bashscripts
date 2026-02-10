#!/bin/bash -l

#SBATCH --time=2:00:00
#SBATCH --account=def-cdeboer
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --output=logs/4b_concatanate_fastas_%j.out
#SBATCH --error=logs/4b_concatanate_fastas_%j.err

HUMAN_FA=data/hg38.fa
YEAST_FA=data/sacCer3.fa
# This uses sed to find the '>' at the start of every header 
# and replaces it with '>yeast_'
sed 's/>/>yeast_/g' ${YEAST_FA} > "$YEAST_FA fixed"
(cat "$YEAST_FA fixed"; echo ""; cat ${HUMAN_FA}) > data/yeast_human_hybrid.fa