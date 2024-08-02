#!/usr/bin/env python3

import pandas as pd
import numpy as np


def _df_shannon_entropy(counts):
    '''Calculate shannon entropy from an array of counts.

    In

    Parameters
    ----------
    counts : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    '''
    probabilities = counts / counts.sum()
    return -np.sum(probabilities * np.log2(probabilities))


def calculate_shannon_entropy(df: pd.DataFrame,
                              group_col: str,
                              value_col: str,
                              result_col: str ='shannon_entropy') -> pd.DataFrame:
    '''Calculate Shannon entropy of count distribution for each group in a pandas DataFrame.

    Assumes that counts have already been computed for each unique value in the group

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    group_col : str
        column name for grouping
    value_col : str
        olumn name for count values to calculate entropy
    result_col : str, optional
        column name for storing calculated entropy, by default 'shannon_entropy'

    Returns
    -------
    pd.DataFrame
        pandas DataFrame with group values and corresponding Shannon entropy
    '''
    result = df.groupby(group_col)[value_col].agg(_df_shannon_entropy).reset_index(name=result_col)

    return result

# # Example usage
# data = {'id_col': [1, 1, 1, 2, 2, 3, 3],
#         'overhang_len': [2, 3, 4, 1, 2, 2, 3],
#         'read_count': [10, 5, 5, 8, 12, 10, 10]}

# # Expected output:
# #    id_col  shannon_entropy
# # 0       1         1.500000
# # 1       2         0.970951
# # 2       3         1.000000

# df = pd.DataFrame(data)

# result_df = calculate_shannon_entropy(df, 'id_col', 'read_count')

# print(result_df)
# #    id_col  shannon_entropy
# # 0       1         1.500000
# # 1       2         0.970951
# # 2       3         1.000000
