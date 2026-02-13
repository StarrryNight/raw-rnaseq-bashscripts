#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --account=rrg-cdeboer
#SBATCH --job-name=bamcoverage
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --output=logs/5_create_bw_%j.out
#SBATCH --error=logs/5_create_bw_%j.err



# Update these to match your actual directory structure
SAMPLE="yac_sample_1_"
BAM="outputs/alignment_2/${SAMPLE}Aligned.sortedByCoord.out.bam"
OUTDIR="outputs/bigwig/${SAMPLE}"
mkdir -p ${OUTDIR}

module load samtools
samtools index ${BAM}
# Generate Forward Strand BigWig
bamCoverage -b ${BAM} \
    -o ${OUTDIR}/${SAMPLE}_forward.bw \
    -of bigwig \
    --filterRNAstrand forward \
    --binSize 1 \
    --normalizeUsing CPM \
    --numberOfProcessors 8

# Generate Reverse Strand BigWig
bamCoverage -b ${BAM} \
    -o ${OUTDIR}/${SAMPLE}_reverse.bw \
    -of bigwig \
    --filterRNAstrand reverse \
    --binSize 1 \
    --normalizeUsing CPM \
    --numberOfProcessors 8