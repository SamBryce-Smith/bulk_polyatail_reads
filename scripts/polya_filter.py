#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
from features import calculate_shannon_entropy
from polya_cluster import assign_id


def parquet_to_gr(parquet_path: str,
                        row_filters: list,
                        return_rle: bool = True):
    '''Read in a parquet file, applying row level filters and returning a pr.PyRanges / pyrle object

    Parameters
    ----------
    parquet_path : str
        _description_
    read_level_filters : list
        _description_
    return_rle : bool, optional
        _description_, by default True
    '''
    gr = pr.PyRanges(pd.read_parquet(parquet_path, filters=row_filters, engine="pyarrow"))

    if return_rle:
        return gr.to_rle()
    else:
        return gr

def vlasenok_filter(parquet_path: str | list,
                    read_level_filters: list = [('overhang_len', '>=', 6), ('frac_As', '>=', 0.8)],
                    shannon_entropy_cutoff: float = 2.0):
    '''Implement Vlasenok et al. filtering criteria

    Parameters
    ----------
    parquet_path : str
        _description_
    read_level_filters : list, optional
        _description_, by default [('overhang_len', '>=', 6), ('frac_As', '>=', 0.8)]
    shannon_entropy_cutoff : float, optional
        _description_, by default 2.0

    Returns
    -------
    _type_
        pyrle object of genome-wide coverage (pooled if multiple samples)
    '''

    if isinstance(parquet_path, str):
        print("Reading from parquet file, applying read-level filters...")
        gr = pr.PyRanges(pd.read_parquet(parquet_path, filters=read_level_filters, engine="pyarrow"))
        
        # assign identifier column (position_id) for each position
        gr = assign_id(gr)
    
    else:
        # list of paths - need to filter all samples, then pool across datasets
        i = 1
        num_samples = len(parquet_path)
        grs = []
        for p in parquet_path:
            print(f"Filtering at read-level for file {i} / {num_samples}")
            tmp_gr = parquet_to_gr(p, read_level_filters, return_rle=False)
            print(f"number of reads - {len(tmp_gr)}")
            # assign identifier column (position_id) for each position
            tmp_gr = assign_id(tmp_gr)
            grs.append(tmp_gr)
            i +=1

        # combine across samples
        # TODO: For large sample sizes this is likely to be very memory intensive - should consider doing the minimal collapsing for shannon entropy early
        gr = pr.concat(grs)
        print(f"Pooled number of reads - {len(gr)}")


    # Calculate shannon entropy per position
    print("Calculating shannon entropy...")
    n_ids = set(gr.as_df()["position_id"])
    print(f"Number of positions for which calculating shannon entropy - {len(n_ids)}")

    # current gr has each individual read per row
    # shannon entropy fnc requires df of id | overhang_len | count where count = number of reads with that overhang_len
    # TODO: try gr.as_df()[["position_id", "overhang_len"]].value_counts() ? Might be more efficient?
    shannon_df = gr.as_df().groupby(["position_id", "overhang_len"]).size().reset_index(name="read_count")

    # calculate shannon entropy for each position (adds shannon_entropy column)
    shannon_df = calculate_shannon_entropy(shannon_df, group_col="position_id", value_col="read_count")

    # filter for IDs passing threshold
    shannon_pass_ids = set(shannon_df.loc[shannon_df["shannon_entropy"].ge(shannon_entropy_cutoff), "position_id"])
    print(f"Number of positions passing overhang shannon entropy threshold of {shannon_entropy_cutoff} - {len(shannon_pass_ids)}")

    # subset gr for positions passing entropy cutoff
    gr = gr.subset(lambda df: df["position_id"].isin(shannon_pass_ids))

    return gr.to_rle()


def twoclass_length_filter(parquet_path: str,
                           read_level_filters: list = [[('overhang_len', '<', 6), ('overhang_len', '>=', 3), ('frac_As', '==', 1)], [('overhang_len', '>=', 6), ('frac_As', '>=', 0.8)]]):
    '''Apply two-class filter based on overhang length.
    
    Two possible ways for events to pass filter:
    1. Below a certain length = (strict) A content threshold
    2. Above a certain length = (relaxed) A content threshold

    Default is a simple two stage filter (less than 6 & 100 % A OR >= 6 & >= 80 % As.

    Parameters
    ----------
    parquet_path : str
        _description_
    read_level_filters : list, optional
        _description_, by default [('overhang_len', '>=', 6), ('frac_As', '>=', 0.8)]
    
    Returns
    -------
    _type_
        pyrle object of genome-wide coverage
    '''

    # read in junction reads from parquet file, applying desired filters. Then convert to rle to represent coverage
    return parquet_to_gr(parquet_path, read_level_filters, return_rle=True)