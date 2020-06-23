import re, math
import pandas as pd
import numpy as np


def read(filename):
    """
    Read an XPS file into a pandas dataframe
    :param filename:
    :return:
    """
    names, header_row = process_header(filename)
    print("Column names found in the {}".format(filename), names)


    res = pd.read_csv(filename, delim_whitespace=True, skip_blank_lines=True,
                       names=names, skiprows=header_row+1)
    print("Read from {}\n".format(filename), res)
    return res


def process_header(filename):
    """
    Hard-coded section to process the header region of an XPS file.
    1. Copy the header line as variable names
    2. Remove all redundant information in the header line via regex substitution

    Returns list of column names, and the zero-based line_number
    where the header line is located

    This particular program is tailored for Jin_O1s_original/100mtorr-150degree.txt,
    one of the results from Ag experiments.
    """
    line_num = 3

    # string = "KE_O1s 670eV	BE_O1s 670eV	CPS_O1s 670eV	O1s 670eV_1_O1s 670eV	O1s 670eV_2_O1s 670eV	O1s 670eV_3_O1s 670eV	O1s 670eV_4_O1s 670eV	O1s 670eV_5_O1s 670eV	Background_O1s 670eV	Envelope_O1s 670eV"
    string = ""
    with open(filename) as fp:
        for i, line in enumerate(fp):
            if i == line_num:
                string = line
    new_str = re.sub('\s670eV_*', '', string)
    new_str = re.sub('[_\s]O1s', ' ', new_str)
    return new_str.split(), line_num


def int_col_names(df):
    """
    Return a new column where all the digit-based column names are changed from str type
    to int type.
    """
    mapper = lambda n: int(n) if isinstance(n, str) and str.isdigit(n) else n
    return df.rename(mapper, axis="columns")


def subtract_background(df):
    """
    Subtract the "Background" column from all columns named by a numerical value or named "Envelope"

    Drop the "Background" column
    """
    res = pd.DataFrame()
    res["BE"] = df["BE"]
    for col in df.columns:
        if isinstance(col, int) or col == "Envelope":
            res[col] = df[col] - df["Background"]
    return res


def read_and_process(filename, round=-1, be_as_index=True):
    """
    Read from filename and process the read dataframe.
    """
    res = read(filename)

    res = subtract_background(int_col_names(res))

    if round > 0:
        res["BE"] = res["BE"].round(round)
        res = res.groupby(['BE']).mean().reset_index()
    if be_as_index:
        res = res.set_index(["BE"])

    # print("Processed from {}".format(filename))
    # with pd.option_context('display.max_rows', 1000, 'display.max_columns', None):
    #     print(res)
    return res


def process_and_pickle(filename):
    """
    Read, process, and pickle the resulting dataframe
    """
    read_and_process(filename).to_pickle("{}_pickled".format(filename))


# process_and_pickle("Jin_O1s_original/100mtorr-150degree.txt")


def find_largest_and_peaks(*files):
    max_val, best_col, best_be, best_file = -math.inf, "", 0, ""
    for f in files:
        print("Reading file {}".format(f))

        df = read_and_process(f)

        for col in df.columns:
            if col != "Envelope" and col != "BE":

                print(col)

                col_max = df.loc[df[col].idxmax()]
                if col_max[col] > max_val:
                    max_val = col_max[col]
                    best_col = col
                    best_be = col_max["BE"]
                    best_file = f
    return max_val, best_col, best_be, best_file

# print(find_largest_and_peaks("Jin_O1s_original/10-2 mtorr-25degree.txt",
#                         "Jin_O1s_original/30 mtorr-25degree.txt",
#                         "Jin_O1s_original/100mtorr-300degree.txt"))

#### Utility methods ####
def fill_zeros(df, start= None, end=None, step=0.1):
    if end:
        r = np.arange(df.index.max() + step, end + step, step)

        zeros = pd.DataFrame(0, index=r, columns=["Intensity"])
        zeros = zeros.reindex(index=zeros.index[::-1])

        df = pd.concat([zeros, df], axis=0)
        # df["BE"] = df.index
    if start:
        r = np.arange(start, df.index.min(), step)

        zeros = pd.DataFrame(0, index=r, columns=["Intensity"])
        zeros = zeros.reindex(index=zeros.index[::-1])

        df = pd.concat([df, zeros], axis=0)
        # df["BE"] = df.index
    return df
