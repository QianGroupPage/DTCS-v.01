import re

from lblcrn.xps_data_processing.condition import Condition
from lblcrn.xps_data_processing.utilities import read_line_blocks


def read_digital_notebook(filename):
    """
    Every line block before the first line block containing a measurement entry will be the general comments section.

    Afterwards, every line block not containing a measurement will be skipped. Each line block will be considered a
    measurement condition.

    :param filename: path to the digital notebook.
    :return: a list of general comments line blocks, and a list of conditions described by the notebook.
    """
    condition_line_blocks = read_line_blocks(filename)
    if condition_line_blocks[0][0].startswith("\ufeffa"):
        condition_line_blocks[0][0] = condition_line_blocks[0][0][1:]
        print(f"Ignoring '\\ufeff' at beginning of file {filename}. " +
              f"This redundant character was attached by the Windows NotePad.")

    general_comments_block_ending = 0
    for i, condition_line_block in enumerate(condition_line_blocks):
        if is_measurement_entry(condition_line_block):
            general_comments_block_ending = i
            break

    general_comments_line_blocks = condition_line_blocks[:general_comments_block_ending]

    conditions = []
    data_line_blocks = condition_line_blocks[general_comments_block_ending:]
    for line_block in data_line_blocks:
        # Skip comment line-blocks
        if not is_measurement_entry(line_block):
            print("skipping comment entry")
            continue
        condition = Condition(line_block)
        conditions.append(condition)
    return general_comments_line_blocks, conditions


def is_measurement_entry(line_block):
    """
    :param line_block: a list of lines, each line represented as a string.
    :return: True if the line-block contains an entry for measurement.
    """
    for line in line_block:
        # Select the lines representing one measurement.
        if re.match('[0-9]+.*', line):
            return True
