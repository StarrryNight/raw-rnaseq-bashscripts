import pyBigWig
import numpy as np
import argparse

def bigwig_to_numpy(bw_path, chrom, start, end, bin_size=1):
    # Open the BigWig
    bw = pyBigWig.open(bw_path)
    
    # Extract values for the range
    # 'numpy=True' makes it significantly faster
    values = bw.values(chrom, start, end, numpy=True)
    
    # Handle NaNs (BigWig uses NaN for zero-coverage regions)
    values = np.nan_to_num(values)
    
    # If binning is required (e.g., 128bp bins for Yorzoi)
    if bin_size > 1:
        # Reshape and average
        n_bins = len(values) // bin_size
        values = values[:n_bins * bin_size].reshape(n_bins, bin_size).mean(axis=1)
    
    bw.close()
    return values

if __name__ == "__main__":
    # Settings based on your recent YAC mapping
    argparser = argparse.ArgumentParser(description='Convert BigWig coverage to NumPy array')
    argparser.add_argument('--chrom', type=str, required=True, help='chromosome name (e.g., chr1)')
    argparser.add_argument('--start', type=str, required=True, help='start of chromosome')
    argparser.add_argument('--end', type=str, required=True, help='end of chromosome')
    argparser.add_argument('--fwd_bw', type=str, required=True)
    argparser.add_argument('--rev_bw', type=str, required=True)
    argparser.add_argument('--output_dir', type=str, required=True, help='output directory path')
    args = argparser.parse_args()
    CHROME = args.chrom
    START = int(args.start)
    END = int(args.end)
    
    # Convert Forward Strand
    fwd_array = bigwig_to_numpy(args.fwd_bw, CHROME, START, END)
    np.save(f"{args.output_dir}/yac_sample_1_fwd.npy", fwd_array)
    
    # Convert Reverse Strand
    rev_array = bigwig_to_numpy(args.rev_bw, CHROME, START, END)
    np.save(f"{args.output_dir}/yac_sample_1_rev.npy", rev_array)
    
    print(f"Success! Saved arrays with shape: {fwd_array.shape}")