import matplotlib.pyplot as plt
import pandas as pd
import argparse


def main():
    argparser = argparse.ArgumentParser(description='Plot YAC coverage data')
    argparser.add_argument('--input', type=str,required=True, help='')
    argparser.add_argument('--output', type=str, required=True, help='Output image file')
    args = argparser.parse_args()
    df = pd.read_csv(args.input, sep='\t', names=['chr', 'pos', 'depth'])
    plt.plot(df['pos'], df['depth'])
    plt.title('YAC Read Depth on Chromosome 1')
    plt.xlabel('Position (bp)')
    plt.ylabel('Coverage Depth')
    plt.show()
    plt.savefig(args.output, dpi=300)
if __name__ == "__main__":
    main()