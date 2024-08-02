#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
import polya_filter as flt
from functools import reduce
from operator import add
from polya_cluster import assign_id, cluster_polya_sites
import argparse
import sys
import os

'''
Script to filter and cluster polyA junction reads into polyA site clusters
1. Apply read-level filters to define polyA junction containing reads
2. (Optionally pooling across samples), cluster closely spaced reads into polyA site clusters (using PolyASite algorithm)
3. Output a BED file of polyA site coordinates, containing the best supported position, read counts etc.
'''


# Each sample, filter from parquet files & convert to pyranges object, assign ID for position
# run clustering function, output BED

# keys for different default filtering criteria
filter_methods = ['vlasenok', 'two_class_simple', "vlasenok_two_class"]



def main(parquet_paths: list,
         filter_method: str,
         read_level_filters: list,
         entropy_cutoff: float,
         window_len: int,
         cluster_per_sample: bool,
         parquet_suffix: str,
         nb_cpu: int,
         output_prefix: str,
         valid_method_keys: list = filter_methods):
    '''_summary_

    Parameters
    ----------
    parquet_paths : list
        _description_
    filter_method : str
        _description_
    window_len : int
        _description_
    cluster_per_sample : bool
        _description_
    parquet_suffix : str
        _description_
    nb_cpu : int
        _description_
    output_prefix : str
        _description_
    '''

    assert not (filter_method in ["vlasenok", "vlasenok_two_class"] and cluster_per_sample), "Clustering per-sample is incompatible with the vlasenok filter methods"
    assert filter_method in valid_method_keys

    # Some methods require pre-filtering/pooling
    if filter_method in ["vlasenok", "vlasenok_two_class"]:
        polya_jncs = flt.vlasenok_filter(parquet_paths, read_level_filters=read_level_filters,shannon_entropy_cutoff=entropy_cutoff)
        # read_level_filters=[[('overhang_len', '<', 6), ('overhang_len', '>=', 3), ('frac_As', '==', 1)], [('overhang_len', '>=', 6), ('frac_As', '>=', 0.8)]
        
        # filtering method already pools across samples (i.e. ready for clustering)
        pooled_cov = polya_jncs

    # remaining filter per-sample only
    elif filter_method == "two_class_simple":

        # extract sample names from .parquet paths
        sample_names = [os.path.basename(p).removesuffix(parquet_suffix) for p in parquet_paths]

        # dict storing polyA junction read coverage in pyrle objects {sample_name: cov_rle}
        sample2cov = {}
    
        i = 1
        num_samples = len(sample_names)

        for path, sample_name in zip(parquet_paths, sample_names):
            print(f"Importing and filtering sample {sample_name} - {i} / {num_samples}")

            # read in junction reads from parquet file, applying desired filters. Then convert to rle to represent coverage
            if filter_method == "two_class_simple":
                polya_jncs = flt.twoclass_length_filter(path, read_level_filters=read_level_filters)

            if cluster_per_sample:
                print(f"Computing polyA site clusters for sample {sample_name}")
                
                cov = polya_jncs.to_ranges()
                if not isinstance(cov, pr.PyRanges):
                    cov = pr.PyRanges(cov)

                # Add id column for each unique position in the genome covered by a polyA read cluster
                cov = assign_id(cov, out_col="Name")

                clusters = cluster_polya_sites(cov, id_col="Name", count_col="Score", cluster_width=window_len, nb_cpu=nb_cpu)
                
                print("Outputting sample-level polyA clusters to file...")
                clusters.to_bed(".".join([output_prefix, sample_name, "polya_clusters.bed"]))

            sample2cov[sample_name] = polya_jncs
            i +=1

        # print(sample2cov[sample_names[0]])
        # pool coverage across samples
        print("Pooling coverage across replicates...")
        pooled_cov = reduce(add, sample2cov.values())
        # pooled cov is dict of {(chr, strand) : Rle}
        # Rle.values accesses the coverage value for given runs in genome
        # sum coverage for each chrom/strand pair, then sum across all pairs
        print(f"Total number of pooled reads - {sum(np.sum(v.values) for _, v in pooled_cov)}")    
    
    # print(pooled_cov.to_ranges())
    # print(sample2cov[sample_names[0]].to_ranges())

    print("Clustering pooled polyA junction reads into polyA site clusters...")
    pooled_cov = pooled_cov.to_ranges()
    if not isinstance(pooled_cov, pr.PyRanges):
        pooled_cov = pr.PyRanges(pooled_cov)

    # Add id column for each unique position in the genome covered by a polyA read cluster
    pooled_cov = assign_id(pooled_cov, out_col="Name")

    pooled_clusters = cluster_polya_sites(pooled_cov, id_col="Name", count_col="Score", cluster_width=window_len, nb_cpu=nb_cpu)
    # print(pooled_clusters)

    # make sure Score column is int dtype (as it's just counts)
    pooled_clusters.Score = pooled_clusters.Score.astype(np.int64)

    print(f"Writing BED file of pooled polyA clusters to file - {output_prefix + '.bed'}")
    pooled_clusters.to_bed(output_prefix + ".bed")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter and cluster polyA junction reads into a BED file of polyA site clusters',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                    )

    parser.add_argument('-i', '--input-parquet', type=str, nargs="+", default=argparse.SUPPRESS, help='Path to .parquet file(s) of extracted soft-clipped reads (produced by extract_polya_reads.py). If multiple, pass consecutively and space-separated')
    parser.add_argument('-f', '--filter-method', choices=filter_methods, help="Method to use to filter for valid polyA junction reads/positions", default="two_pass_simple")
    parser.add_argument('-w','--cluster-window', type=int, default=12, help='Length of interval to extend either side of provided PAS coordinates to cluster polyA junction reads (note that polyA site coordinates will be 1nt long, so total max length of cluster = 2*window + 1)')
    parser.add_argument('-s','--suffix', type=str, default=".parquet", help='Suffix to strip from input parquet file path(s) to generate sample name')
    parser.add_argument("-e", "--entropy-cutoff", type=float, default=2.0, help="Minimum shannon entropy threshold for valid positions. Only applies to 'vlasenok*' filter methods")
    parser.add_argument("--long-overhang-min-length", type=int, default=6, help="Minimum overhang length to define as 'long' overhangs")
    parser.add_argument("--short-overhang-length-range", type=str, default="3-5", help="Range of overhang lengths to define 'short' overhangs. 3-5 = 3,4,5 (i.e. includes both boundaries of interval). Only applies to '*two_class' filter methods.")
    parser.add_argument("--long-overhang-fraction-a", type=float, default=0.8, help="minimum A/tail nucleotide content for 'long' overhang reads to be considered valid")
    parser.add_argument("--short-overhang-fraction-a", type=float, default=1.0, help="Range of A/tail nucleotide content for 'short' overhang reads to be considered valid. Only applies to '*two_class' filter methods.")
    parser.add_argument("--per-sample", action="store_true", help="Whether to additionally output clusters per sample")
    parser.add_argument("-c", "--cores", type=int, default=1, help="Number of cores to use for parallel processing (currently just for PAS clustering)")
    parser.add_argument('-o','--output_prefix', type=str, default="pas_clusters", help='Prefix for output file names (<output_prefix>.bed)')

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    # construct read level parquet filters from input thresholds
    filt_long_overhang = [('overhang_len', '>=', args.long_overhang_min_length), ('frac_As', '>=', args.long_overhang_fraction_a)]
    
    if args.filter_method in ['two_class_simple', "vlasenok_two_class"]:
        short_len_range = args.short_overhang_length_range.split("-") # 3,5 -> [3,5]
        assert len(short_len_range) == 2

        short_overhang = [[('overhang_len', '<=', int(short_len_range[1])), ('overhang_len', '>=', int(short_len_range[0])), ('frac_As', '>=', args.short_overhang_fraction_a)]]

        # append to keep filter tuples as a list
        short_overhang.append(filt_long_overhang)

        read_level_filters = short_overhang

    else:
        # using a single read length + a content filter
        read_level_filters = filt_long_overhang


    print(read_level_filters)


    main(args.input_parquet, args.filter_method, read_level_filters, args.entropy_cutoff, args.cluster_window, args.per_sample, args.suffix, args.cores, args.output_prefix)
