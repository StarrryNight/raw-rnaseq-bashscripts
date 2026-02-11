

module load samtools
HUMAN_FASTA="data/hg38.fa"
YEAST_FASTA="data/sacCer3.fa"


read TOP_CHR START END < my_yac.bed

# 1. Extract the human slice (using the coordinates you found)
HUMAN_YAC_INSERT=$(samtools faidx $HUMAN_FASTA $TOP_CHR:$START-$END) 


# 3. Combine with the Yeast genome
cat $YEAST_FASTA $HUMAN_YAC_INSERT > yeast_human_chimeric.fa