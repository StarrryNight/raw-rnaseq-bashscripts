#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --account=rrg-cdeboer
#SBATCH --job-name=bamcoverage
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --output=logs/7_generate_numpy_%j.out
#SBATCH --error=logs/7_generate_numpy_%j.err

read CHROME START END < my_yac.bed
echo "Generating NumPy array for ${CHROME}:${START}-${END}"
python bw_to_np.py --chrom=${CHROME} --start=${START} --end=${END} \
    --fwd_bw=outputs/bigwig/yac_sample_1/yac_sample_1_forward.bw \
    --rev_bw=outputs/bigwig/yac_sample_1/yac_sample_1_reverse.bw \
    --output_dir=final_output
