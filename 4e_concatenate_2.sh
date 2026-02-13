#!/bin/bash
module load samtools

# 1. Setup Paths
HUMAN_FASTA="data/hg38.fa"
YEAST_FASTA="data/sacCer3.fa"
BED_FILE="my_yac.bed"
OUTPUT_FASTA="yeast_human_chimeric.fa"

# 2. Extract variables using awk (much more reliable than 'read')
# This handles tabs or spaces automatically
TOP_CHR=$(awk '{print $1}' "$BED_FILE")
START=$(awk '{print $2}' "$BED_FILE")
END=$(awk '{print $3}' "$BED_FILE")

echo "Targeting Region: $TOP_CHR from $START to $END"

# 3. Validation: If TOP_CHR is empty, stop now!
if [ -z "$TOP_CHR" ]; then
    echo "ERROR: Could not read my_yac.bed. Is the file empty?"
    exit 1
fi

# 4. Extract and build
# We use a subshell to rename the header on the fly to avoid temp files
samtools faidx "$HUMAN_FASTA" "${TOP_CHR}:${START}-${END}" | \
sed "s/^>.*/>${TOP_CHR}/" > human_yac_insert.fa

# Check if the extraction actually worked
if [ ! -s human_yac_insert.fa ]; then
    echo "ERROR: Human extraction failed. Check if $TOP_CHR exists in $HUMAN_FASTA"
    exit 1
fi

# 5. Combine
cat "$YEAST_FASTA" human_yac_insert.fa > "$OUTPUT_FASTA"

# 6. Index the new genome
samtools faidx "$OUTPUT_FASTA"

echo "Success! Final file size: $(du -sh $OUTPUT_FASTA)"