#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np

# test_dict = {"Chromosome": ["chr1"]*6,
#  "Start": [50,55,60,65,80,88],
#  "End": [51,56,61,66,81,89],
#  "Strand": ["+"]*6,
#  "read_count": [5,1,1,1,1,2]}

# test_gr = pr.from_dict(test_dict)

# # print(test_gr)

# # Cluster with a slack of 25, are they all grouped together?
# test_gr = test_gr.cluster(strand=True, slack=25)
# print(test_gr)

# Despite 88 being > 25 from 50 they are grouped together.
# If want a maximum cluster size, need additional filtering

# PolyASite clustering approach
# Extend every read by 24nt (so 25nt window)
# Sort by highest read support, collapse more weakly expressed sites to the highest supported event provided...
# within 24 nt window around higher expressed site)
# Iterate until all reads assigned to cluster

# sort dfs by cluster and decreasing read count
# test_gr = test_gr.apply(lambda df: df.sort_values(by=["Cluster", "read_count"], ascending=False))
# # print(test_gr)

# # Need a position ID (chr, start, end, strand)
# # function applied to each df, 
# test_gr = test_gr.assign("position_id",
#                          lambda df: df["Chromosome"].str.cat(df[["Strand", "Start", "End"]].astype(str), sep=":"))

# print(test_gr)

# Use pr.cluster as an initial clustering (efficient but imprecise)
# collapse within cluster
# repeat until complete 

def assign_id(gr: pr.PyRanges,
              out_col: str = "position_id",
              cols_to_cat: list = ["Chromosome", "Strand", "Start", "End"],
              sep: str = ":"):
    '''Assign a unique identifier for each interval

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    out_col : str, optional
        _description_, by default "position_id"
    cols_to_cat : list, optional
        _description_, by default ["Chromosome", "Strand", "Start", "End"]
    sep : _type_, optional
        _description_, by default ":"
    '''

    assert all(col in gr.columns for col in cols_to_cat)

    return gr.assign(out_col,
                     lambda df: df[cols_to_cat[0]].str.cat(df[cols_to_cat[1:]].astype(str), sep=sep))



def select_first_value(series: pd.Series):
    return series.iat[0]



def _df_collapse_cluster(df: pd.DataFrame,
                     id_col: str = "position_id",
                     count_col: str = "read_count",
                     max_distance: int = 25
                     ):
    
    # find the position with strongest read support within cluster (representative position)
    # Assuming position sorted, most proximal position arbitrarily selected in case of ties 
    best_idx = df[count_col].idxmax()

    # Extract values for Chromosome, Strand, and id_col based on best_idx
    chromosome = df.at[best_idx, "Chromosome"]
    strand = df.at[best_idx, "Strand"]
    id_value = df.at[best_idx, id_col]
    best_start = df.at[best_idx, "Start"]

    # Filter for reads falling within max_distance nt of best supported position
    df_filt = df[np.abs(df["Start"].values - best_start) <= max_distance]

    # Define output columns for collapsed cluster entry
    agg_result = {
        "Chromosome": chromosome,
        "Strand": strand,
        "Start": df_filt["Start"].min(),
        "End": df_filt["End"].max(),
        count_col: df_filt[count_col].sum(),
        id_col: id_value
    }
    
    # Convert aggregation result to DataFrame
    df_out = pd.DataFrame([agg_result])

    # Convert coordinate columns to integers
    df_out = df_out.astype({"Start": int, "End": int})

    return df_out

# print(_df_collapse_cluster(test_gr.as_df()))


def cluster_polya_sites(gr: pr.PyRanges,
                        id_col: str = "position_id",
                        count_col: str = "read_count",
                        cluster_width: int = 25,
                        nb_cpu: int = 1):
    '''_summary_

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    id_col : str, optional
        _description_, by default "position_id"
    count_col : str, optional
        _description_, by default "read_count"
    cluster_width : int, optional
        _description_, by default 25
    '''

    gr_cluster = gr.cluster(strand=True, slack=cluster_width)
    gr_cluster_collapsed = pr.PyRanges()
    i = 0
    # Workflow = rounds of collapsing clusters until no more unassigned reads
    #1. first round cluster + collapse
    #2. negated overlap initial gr with existing clusters,

    while True:
        if i == 0:
            clpsd = gr_cluster.apply(lambda df: df.groupby("Cluster").apply(lambda df: _df_collapse_cluster(df, id_col=id_col, count_col=count_col, max_distance=cluster_width)), nb_cpu=nb_cpu)
            gr_cluster_collapsed = clpsd
            i+=1 
            # now try again
            continue

        # check for overlaps between original gr and collapsed clusters
        gr_to_clps = gr.overlap(gr_cluster_collapsed, strandedness="same", invert=True)

        if len(gr_to_clps) == 0:
            # all reads assigned to clusters, stop collapsing
            break

        # another round of collapsing
        gr_to_clps = gr_to_clps.cluster(strand=True, slack=cluster_width)
        clpsd = gr_to_clps.apply(lambda df: df.groupby("Cluster").apply(lambda df: _df_collapse_cluster(df, id_col=id_col, count_col=count_col, max_distance=cluster_width)), nb_cpu=nb_cpu)
        
        # update collapsed clusters
        gr_cluster_collapsed = pr.concat([gr_cluster_collapsed, clpsd])
        i+=1 


    print(f"PolyA site clustering complete - number of iterations required: {i}")
    return gr_cluster_collapsed


# test_dict = {"Chromosome": ["chr1"]*6,
#  "Start": [50,55,60,65,80,88],
#  "End": [51,56,61,66,81,89],
#  "Strand": ["+"]*6,
#  "read_count": [5,1,1,1,1,2]}

# test_gr = pr.from_dict(test_dict)
# test_gr = assign_id(test_gr)
# # print(test_gr)

# # Cluster with a slack of 25, are they all grouped together?
# test_gr = test_gr.cluster(strand=True, slack=25)
# print(test_gr)

# # Despite 88 being > 25 from 50 they are grouped together.
# # If want a maximum cluster size, need additional filtering

# print(cluster_polya_sites(test_gr.drop("Cluster")))