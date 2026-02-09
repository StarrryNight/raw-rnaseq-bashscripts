#!/bin/bash -l

#SBATCH --time=2:00:00
#SBATCH --account=def-cdeboer
#SBATCH --job-name=fastqc_raw
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --output=fastqc_raw_%j.out #SBATCH --error=fastqc_raw_%j.err


module load trimmomatic

#!/bin/bash

# Define paths to your input files
IN1="data/1951f6-DL-8_10_S13_L007_R1_001.fastq.gz"
IN2="data/1951f6-DL-8_10_S13_L007_R2_001.fastq.gz"

# Define output names (I added 'trimmed_' to distinguish them)
OUT1_P="outputs/fastqc/after_trim/trimmed_R1_paired.fq.gz"
OUT1_U="outputs/fastqc/after_trim/trimmed_R1_unpaired.fq.gz"
OUT2_P="outputs/fastqc/after_trim/trimmed_R2_paired.fq.gz"
OUT2_U="outputs/fastqc/after_trim/trimmed_R2_unpaired.fq.gz"

# The command
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
    -threads 4 \
    -phred33 \
    "$IN1" "$IN2" \
    "$OUT1_P" "$OUT1_U" \
    "$OUT2_P" "$OUT2_U" \
    ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
