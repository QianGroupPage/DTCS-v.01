"""
Common Input and Output functionalities to comprehend tables of all formats.

Credits:
Jin Qian, Ye Wang, Andrew Bogdan, Rithvik Panchapakesan, Ethan Crumlin
"""

import pandas as pd


def read_line_delimited_excel(io):
    """
    Read an Excel sheet delimited by empty lines into several Excel sheets.

    Each df has a header at the first row.

    :param io: path to the Excel sheet or a file-like object;
    :return: a list of Dataframes by sequence.
    """
    df = pd.read_excel(io=io)
    empty_regions = []
    prev_is_empty = False
    for index, row in df.iterrows():
        if all([pd.isnull(e) for e in row]):
            if not empty_regions or prev_is_empty is False:
                empty_regions.append([index])
            elif len(empty_regions[-1]) == 1:
                empty_regions[-1].append(index)
            elif index == empty_regions[-1][1] + 1:
                empty_regions[-1][1] += 1
            else:
                empty_regions.append([index])
    # Redundant code to simplify logic to handle no empty region.
    if not empty_regions:
        return [df]
    dfs = []
    prev_j = -1
    for r in empty_regions:
        if len(r) == 1:
            i, j = r[0], r[0]
        elif len(r) == 2:
            i, j = r
        # Add a df if prev_j and i are not both at 0.
        if prev_j + 1 != i:
            new_df = df.iloc[prev_j + 1:i, :].copy()
            new_df.dropna(axis='columns', how='all', inplace=True)
            if dfs:
                new_df.columns = new_df.iloc[0, :]
                new_df = new_df[1:]
            dfs.append(new_df)
        prev_j = j
    if not prev_j + 1 == len(df):
        new_df = df.iloc[prev_j + 1:, :].copy()
        new_df.dropna(axis='columns', how='all', inplace=True)
        if not new_df.empty and dfs:
            new_df.columns = new_df.iloc[0, :]
            new_df = new_df[1:]
        dfs.append(new_df)
    return dfs
