import re

import numpy as np
import pandas as pd

def read(filename):
    """
    Read an XPS file into a pandas dataframe
    :param filename:
    :return:
    """
    names, header_row = process_header()
    print("Column names found in the {}".format(filename), names)


    res = pd.read_csv(filename, delim_whitespace=True, skip_blank_lines=True,
                       names=names, skiprows=header_row+1)
    print("Read from {}\n".format(filename), res)
    return res

def process_header():
    """
    Hard-coded section to process the header region of an XPS file.
    1. Copy the header line as variable names
    2. Remove all redundant information in the header line via regex substitution

    Returns list of column names, and the zero-based line_number
    where the header line is located

    This particular program is tailored for Jin_O1s/100mtorr-150degree.txt,
    one of the results from Ag experiments.
    """
    line_num = 3
    string = "KE_O1s 670eV	BE_O1s 670eV	CPS_O1s 670eV	O1s 670eV_1_O1s 670eV	O1s 670eV_2_O1s 670eV	O1s 670eV_3_O1s 670eV	O1s 670eV_4_O1s 670eV	O1s 670eV_5_O1s 670eV	Background_O1s 670eV	Envelope_O1s 670eV"
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

def read_and_process(filename):
    """
    Read from filename and process the read dataframe.
    """
    res = read(filename)

    res = subtract_background(int_col_names(res))

    print("Processed from {}".format(filename))
    with pd.option_context('display.max_rows', 1000, 'display.max_columns', None):
        print(res)
    return res

def process_and_pickle(filename):
    """
    Read, process, and pickle the resulting dataframe
    """
    read_and_process(filename).to_pickle("{}_pickled".format(filename))


# process_and_pickle("Jin_O1s/100mtorr-150degree.txt")

first_line = '[region 2]'
start_info = '[info'
start_data = '[data'
info_keys = {'region name', 'center energy'}


def read_data(filename):
    f = open(filename)

    it = file_iter(f)
    line = next(it, None)

    data = []
    # Skip the survey scan
    while line is not None and line.strip().lower() != first_line:
        line = next(it, None)

    while line is not None:
        data.append(parse_region(it))
        line = next(it, None)

    f.close()
    return data

def read_exp(filename, be=0):
    f = open(filename)

    it = file_iter(f)
    next(it, None)
    bes, intensities = parse_new_data(it, be)
    bes.reverse()
    intensities.reverse()
    series = pd.Series(data=intensities, index=bes)
    return series

def read_new_data(filename, be=0):
    f = open(filename)

    it = file_iter(f)
    next(it, None)
    x, y = parse_new_data(it, be)
    x.reverse()
    y.reverse()
    return [XPSData({}, x, y)]


def file_iter(f):
    yield from f


def parse_region(it):
    info = {}
    x, y = [], []
    line = next(it, None)

    # Skip to info section
    while line is not None:
        line = line.strip().lower()
        if start_info in line:
            break

        line = next(it, None)

    while line is not None:
        line = line.strip().lower()
        # Parse a region
        if start_data in line:
            x, y = parse_data(it)
            break

        if line[:line.find('=')] in info_keys:
            info[line[:line.find('=')]] = line[line.find('=')+1:]

        line = next(it, None)

    return XPSData(info, x, y)


def parse_data(it):
    x, y = [], []
    line = next(it, None)

    while line is not None:
        line = line.strip()
        if not line:
           break

        nums = line.split()

        x.append(float(nums[0]))
        y.append(float(nums[1]))
        line = next(it, None)

    return x, y


def parse_new_data(it, be):
    x, y = [], []
    line = next(it, None)

    while line is not None:
        line = line.strip()
        if not line:
           break

        nums = line.split()

        x.append(float(nums[be]))
        y.append(float(nums[8]) - float(nums[7]))
        line = next(it, None)

    return x, y


class XPSData:
    def __init__(self, info, binding_energies, intensities):
        self.info = info
        self.binding_energies = np.asarray(binding_energies)
        self.intensities = np.asarray(intensities)
