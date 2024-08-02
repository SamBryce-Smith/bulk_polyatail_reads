#!/usr/bin/env python3

import argparse
import pandas as pd
from functools import reduce
import sys

# merge two DataFrames on the common ID column
def merge_dfs(df1, file_path, id_column):
    
    df2 = pd.read_csv(file_path, sep='\t')
    
    if len(df1) == 0:
        return df2
    else:  
        return pd.merge(df1, df2, on=id_column)


def main():

    parser = argparse.ArgumentParser(description='Merge multiple polyA junction count matrices (produced by count_pas.py) into a single matrix TSV files into single matrix')

    # Add arguments
    parser.add_argument('input_files', nargs='+', help='List of input TSV files to merge')
    parser.add_argument('-o', '--output_file', metavar='output_file', default='merged_output.tsv',
                        help='Output TSV file name (default: merged_output.tsv)')
    parser.add_argument('-id', '--id_column', metavar='id_column', default='position_id_name',
                        help='Name of the ID column to merge on (default: position_id_name)')
    
    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    # Parse the arguments
    args = parser.parse_args()

    # Use reduce to iteratively merge all DataFrames
    merged_df = reduce(lambda df1, file_path: merge_dfs(df1, file_path, args.id_column), args.input_files, pd.DataFrame())

    # Write the merged DataFrame to the output TSV file
    merged_df.to_csv(args.output_file, sep='\t', index=False)
    print(f'Merged data saved to {args.output_file}')

if __name__ == '__main__':
    main()
