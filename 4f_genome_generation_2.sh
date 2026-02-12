#!/bin/bash -l

#SBATCH --time=2:00:00
#SBATCH --account=rrg-cdeboer
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --output=logs/4f_genome_generation_2_%j.out
#SBATCH --error=logs/4f_genome_generation_2_%j.err

module load star/2.7.11b

GENOME_DIR=data/hybrid_star_index_2
YEAST_HUMAN_FILE=yeast_human_chimeric.fa

mkdir -p ${GENOME_DIR}

STAR --runMode genomeGenerate \
     --runThreadN 8 \
     --genomeDir ${GENOME_DIR} \
     --genomeFastaFiles ${YEAST_HUMAN_FILE} \
     --genomeSAindexNbases 14