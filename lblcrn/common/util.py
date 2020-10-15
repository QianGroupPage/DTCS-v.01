import re

import numpy as np
import pandas as pd
from matplotlib import colors


def weave(*threads):
    n = len(threads)
    woven = np.zeros(np.size(threads))
    for i, t in enumerate(threads):
        woven[i::n] = t
    return woven


def unweave(woven, n):
    return [woven[i::n] for i in range(n)]


def multiple_replace(dict, text):
    """
    :param dict: a dictionary from the string to be replaced to what it should be replaced to
    :return: text with correct replacements.
    """
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))
    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text)


# TODO: I should probably hide these
# TODO: This won't handle colors given by tuples
def scale_color_brightness(color: str, by: float = 1.0) -> str:
    """Scale the given color's brightness by an amount."""
    red, blue, green, alpha = colors.to_rgba(color)
    hue, sat, value = colors.rgb_to_hsv((red, blue, green))
    value *= by
    # Has to be in range [0, 1]
    value = min(value, 1)
    value = max(value, 0)
    red, blue, green = colors.hsv_to_rgb((hue, sat, value))

    return colors.to_hex((red, blue, green, alpha), keep_alpha=True)


# TODO: This won't handle colors given by tuples
def set_color_alpha(color: str, alpha: float = 1.0):
    """Give the given color an alpha value."""
    red, blue, green = colors.to_rgb(color)
    if not 0 < alpha < 1:
        raise ValueError('Alpha must be in range [0, 1]')
    return colors.to_hex((red, blue, green, alpha), keep_alpha=True)


# Dataframe Manipulation
def resample_by_skipping(df, step=1000):
    """
    Resample a dataframe by skipping every thousand or "step" number of rows.
    """
    new_df = pd.DataFrame()
    for i in range(0, len(df.index), step):
        new_df.append(df.iloc[i, :])
    return new_df


