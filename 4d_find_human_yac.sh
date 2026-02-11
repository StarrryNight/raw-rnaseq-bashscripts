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
TOP_CHR=$(samtools idxstats ${YAC_SAMPLE_FILE} | \
grep "^chr" | \
awk '$2 > 10000000' | \
sort -k3,3rn | \
head -n 1 |\
cut -f1)

echo "Top chromosome is: $TOP_CHR"


NUM_READS=$(samtools view -c ${YAC_SAMPLE_FILE} "$TOP_CHR")
START_INDEX=$(( NUM_READS / 100 ))
END_INDEX=$(( NUM_READS * 99 / 100 ))


#TODO rename duplicated chrom names
# This is exactly what you wrote, and it works great 
# because Bash turns the two lines into "number1 number2"
read START END <<< $(samtools view ${YAC_SAMPLE_FILE} "$TOP_CHR" 2>/dev/null | \
    awk '{print $4}' | \
    sort -n | \
    sed -n "${START_INDEX}p;${END_INDEX}p"| \
    tr "\n" " ")

echo "Calculated YAC region: $START to $END"

# 3. Get Depth for the SPECIFIC region
# Instead of sed, we use the -r flag in samtools depth to target the YAC.
# We also downsample (NR % 100) so the Python plot isn't sluggish.
echo "Extracting coverage data..."
samtools depth -a -r "${TOP_CHR}:${START}-${END}" ${YAC_SAMPLE_FILE} | \
    awk 'NR % 100 == 0' > coverage_data.txt


python graph_yac.py --input=coverage_data.txt --output=coverage_plot.png