#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
from features import calculate_shannon_entropy
from polya_cluster import assign_id, cluster_polya_sites
from operator import add
from functools import reduce
import argparse
import sys
import os


'''
Implement the Vlasenok et al. approach for polyA-tail containing reads:
- pool across samples
- overhang >= 6nt long, 80 % A content
- shannon-entropy of overhang length distribution at specific position - cut-off = H >= 2
- cluster remaining positions (for this will just use my implementation of PolyASite)
'''



def main(parquet_paths: list,
         parquet_suffix: str,
         cluster_window_len: int,
         nb_cpu: int,
         output_prefix: str):
    '''_summary_

    Parameters
    ----------
    parquet_paths : list
        _description_
    '''

    # extract sample names from .parquet paths
    sample_names = [os.path.basename(p).removesuffix(parquet_suffix) for p in parquet_paths]

    # dict storing filtered polyA junction read coverage in pyrle objects {sample_name: cov_rle}
    sample2cov = {}
    i = 1
    num_samples = len(sample_names)

    # read-level filters implemented by Vlasenok et al. to extract from parquet file
    vlasenok_readlevel_filt = [('overhang_len', '>=', 6), ('frac_As', '>=', 0.8)]
    shannon_entropy_cutoff = 2

    
    for sample_id, sample_path in zip(sample_names, parquet_paths):

        print(f"Applying read-level filters for sample ID - {sample_id} - ({i} / {num_samples})")
        
        gr = pr.PyRanges(pd.read_parquet(sample_path, filters=vlasenok_readlevel_filt))

        # Calculate shannon entropy per position
        # first assign identifier column (position_id) for each position
        print("Calculating shannon entropy...")
        gr = assign_id(gr)

        n_ids = set(gr.as_df()["position_id"])
        print(f"Number of positions for which calculating shannon entropy - {len(n_ids)}")

        # current gr has each individual read per row
        # shannon entropy fnc requires df of id | overhang_len | count where count = number of reads with that overhang_len
        shannon_df = gr.as_df().groupby(["position_id", "overhang_len"]).size().reset_index(name="read_count")

        # calculate shannon entropy for each position (adds shannon_entropy column)
        shannon_df = calculate_shannon_entropy(shannon_df, group_col="position_id", value_col="read_count")

        # filter for IDs passing threshold
        shannon_pass_ids = set(shannon_df.loc[shannon_df["shannon_entropy"].ge(shannon_entropy_cutoff), "position_id"])
        print(f"Number of positions passing overhang shannon entropy threshold of {shannon_entropy_cutoff} - {len(shannon_pass_ids)}")

        # subset gr for positions passing entropy cutoff
        gr = gr.subset(lambda df: df["position_id"].isin(shannon_pass_ids))

        sample2cov[sample_id] = gr.to_rle()

        i += 1

    print("Filtering for all samples complete")
    print("Pooling coverage across samples...")
    pooled_cov = reduce(add, sample2cov.values())
    # print(pooled_cov.to_ranges())

    print("Clustering pooled polyA junction reads into polyA site clusters...")
    pooled_cov = pooled_cov.to_ranges()
    if not isinstance(pooled_cov, pr.PyRanges):
        pooled_cov = pr.PyRanges(pooled_cov)

    # Add id column for each unique position in the genome covered by a polyA read cluster
    pooled_cov = assign_id(pooled_cov, out_col="Name")

    pooled_clusters = cluster_polya_sites(pooled_cov, id_col="Name", count_col="Score", cluster_width=cluster_window_len, nb_cpu=nb_cpu)
    # print(pooled_clusters)

    print("Writing BED file of pooled polyA clusters to file...")
    pooled_clusters.to_bed(output_prefix + ".all_samples.polya_clusters.bed")

    # TODO: Output stats file with shannon entropy values







        





    