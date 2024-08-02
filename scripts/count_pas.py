#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
from polya_cluster import assign_id
import argparse
import sys
import os

def count_pas(gr: pr.PyRanges,
              parquet_path: str,
              overhang_filters: list = [[('overhang_len', '<', 6), ('frac_As', '==', 1)],
                                        [('overhang_len', '>=', 6), ('frac_As', '>=', 0.8)]],
              filter_chrs: bool = True,
              extend_by: int = 12,
              count_col: str = "read_count"
              ) -> pr.PyRanges:
    '''Assign polyA junction reads to predefined set of polyA sites

    Parameters
    ----------
    gr : pr.PyRanges
        polyAsites to quantify
    parquet_path : str
        _description_
    overhang_filters : list
        how to prefilter junction reads by overhang features. See pd.read_parquet() documentation for how to construct these filters
    filter_chrs : bool, optional
        whether to only load the chromosomes found in gr, by default True
    extend_by : int, optional
        _description_, by default 12
    count_col : str, optional
        Name of output column storing read counts, by default "read_count"
    '''

    # stretch pas by specified distance
    gr_ext = gr.extend(extend_by)

    # read in parquet file to PyRanges object, applying preprovided filters
    # (add chromosome filter if requested)
    if filter_chrs:
        gr_chrs = set(gr.chromosomes)
        chr_filter = [("Chromosome", "in", gr_chrs)]
        # prepend chromosome condition to provided filters
        ohang_filters = [chr_filter + cond for cond in overhang_filters]
    
    else:
        ohang_filters = overhang_filters

    # read in parquet file to pr.PyRanges object
    polya_jncs = pr.PyRanges(pd.read_parquet(parquet_path, filters=ohang_filters))

    # count number of reads for each polyA site/cluster
    gr_ext_counts = gr_ext.count_overlaps(polya_jncs, strandedness="same", keep_nonoverlapping=True, overlap_col=count_col)

    return gr_ext_counts


def count_pas_wrapper(sample_ids: iter, sample_paths: iter, pas_bed: pr.PyRanges, extend_by: int = 24) -> pr.PyRanges:
    '''run the count_pas function across a iterable of paths to generate combined file (with sample of origin information)

    Parameters
    ----------
    sample_ids : iter
        _description_
    sample_paths : iter
        _description_
    pas_bed : pr.PyRanges
        _description_
    extend_by : int, optional
        _description_, by default 24

    Returns
    -------
    pr.PyRanges
        _description_
    '''

    assert len(sample_ids) == len(sample_paths)

    # print(f"Number of samples - {len(sample_ids)}")

    gr_list = []
    i = 1
    num_samples = len(sample_ids)
    
    for sample_id, sample_path in zip(sample_ids, sample_paths):

        print(f"counting for sample ID - {sample_id} - ({i} / {num_samples})")

        # count reads for sample and cryptic PAS
        count_gr = count_pas(pas_bed, sample_path, extend_by=extend_by)

        # assign sample ID as column
        count_gr.sample_name = sample_id

        gr_list.append(count_gr)
        i += 1

    concat_pas_gr = pr.concat(gr_list)

    return concat_pas_gr


def gr_to_count_matrix(gr: pr.PyRanges, region_id_col: str = "position_id_name", sample_id_col: str = "sample_name", count_col: str = "read_count") -> pd.DataFrame:
    '''Generate polyA junction count matrix of interval ids as rows and sample_ids as columns

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    region_id_col : str, optional
        _description_, by default "position_id_name"
    sample_id_col : str, optional
        _description_, by default "sample_name"
    count_col : str, optional
        _description_, by default "read_count"

    Returns
    -------
    pd.DataFrame
        _description_
    '''

    # Convert to pandas dataframe and subset to minimal columns
    df = gr.as_df()[[region_id_col, sample_id_col, count_col]]

    # pivot from long to wide format
    mtx = df.pivot(index=region_id_col, columns=sample_id_col, values=count_col)

    # remove the column index name (sample_id_col, not necessary)
    mtx = mtx.rename_axis(None, axis=1)

    # return the region_id index to a column
    mtx = mtx.reset_index(region_id_col)

    return mtx


def main(parquet_paths: list,
         bed_path: str,
         window_len: int,
         parquet_suffix: str,
         output_prefix: str
         ):
    '''_summary_

    Parameters
    ----------
    parquet_path : list
        _description_
    bed_path : str
        _description_
    window_len : int
        _description_
    parquet_suffix : str
        _description_
    output_prefix : str
        _description_
    '''

    # extract sample names from .parquet paths
    sample_names = [os.path.basename(p).removesuffix(parquet_suffix) for p in parquet_paths]

    # read in and prepare BED
    bed = pr.read_bed(bed_path)
    # assign identifier for each interval (just original coordinates combined)
    bed = assign_id(bed)
    # combine with Name field (can contain useful metadata)
    bed = assign_id(bed, out_col="position_id_name", cols_to_cat=["position_id", "Name"], sep=";")

    # count against provided PAS over input sample names
    counts_long = count_pas_wrapper(sample_ids=sample_names, sample_paths=parquet_paths, pas_bed=bed, extend_by=window_len)

    # generate wide counts matrix (rows = intervals from BED file, columns = sample_names)
    counts_wide = gr_to_count_matrix(counts_long)

    # counts_wide
    print("Outputting wide count matrix to file...")
    counts_wide.to_csv(output_prefix + ".count_matrix.tsv", sep="\t", header=True, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count polyA junction reads overlapping provided BED file of polyA sites',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                      )

    parser.add_argument('-i', '--input-parquet', type=str, nargs="+", default=argparse.SUPPRESS, help='Path to .parquet file(s) of extracted soft-clipped reads (produced by extract_polya_reads.py). If multiple, pass consecutively and space-separated')
    parser.add_argument('-b', '--bed', type=str, default=argparse.SUPPRESS, help='Path to BED file of polyA site coordinates to quantify')
    parser.add_argument('-w','--window', type=int, default=12, help='Length of interval to extend either side of provided PAS coordinates to match polyA junction reads (note that polyA site coordinates should be 1nt long)')
    parser.add_argument('-s','--suffix', type=str, default=".parquet", help='Suffix to strip from input parquet file path(s) to generate sample name')
    parser.add_argument('-o','--output_prefix', type=str, default="pas_counts", help='Prefix for output file names')

    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()
    main(args.input_parquet, args.bed, args.window, args.suffix, args.output_prefix)