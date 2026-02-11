#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --account=rrg-cdeboer
#SBATCH --job-name=bamcoverage
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --output=logs/7_bamcov_%j.out
#SBATCH --error=logs/7_bamcov_%j.err



# Update these to match your actual directory structure
SAMPLE="yac_sample_1"
BAM="outputs/alignment/${SAMPLE}Aligned.sortedByCoord.out.bam"
OUTDIR="outputs/bigwig/${SAMPLE}"
mkdir -p ${OUTDIR}

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