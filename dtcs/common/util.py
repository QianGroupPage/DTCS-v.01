from typing import Union

import functools
import itertools
import importlib
import re
import string
import warnings

import pandas as pd
import numpy as np
import sympy as sym
from matplotlib import colors
from matplotlib import patches as mpatches

from dtcs.common.display import color_map, latex_map


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
        row = df.iloc[i, :]
        row.name = df.index[i]
        new_df.append(row)
    return new_df


def symbol_to_name(sym_or_str):
    """Convert a sym.Symbol or a string to a string through sym.Symbol.name."""
    if type(sym_or_str) == sym.Symbol:
        return sym_or_str.name
    else:
        return sym_or_str


def alphabet():
    """Generate A, B, C, ... Z, AA, AB, ... ZZ, AAA, AAB, and so on."""
    for size in itertools.count(1):
        for s in itertools.product(string.ascii_uppercase, repeat=size):
            yield ''.join(s)


def anon_names(num):
    """Return num anonymous names (from alphabet())."""
    return list(itertools.islice(alphabet(), num))


def flat(lst: Union[list, tuple]):
    """An iterator which flattens a list/tuple of list/tuples."""
    for item in lst:
        if isinstance(item, (list, tuple)):
            yield from flat(item)
        else:
            yield item

def flatten(lst: Union[list, tuple]):
    """Flatten a list/tuple of list/tuples."""
    return list(flat(lst))

# --- Function Decorators -----------------------------------------------------
def depreciate(func):
    @functools.wraps(func)
    def depr_func(*args, **kwargs):
        warnings.warn(f'Depreciated function {func.__module__}.{func.__name__}', DeprecationWarning)
        return func(*args, **kwargs)
    return depr_func


EXTRAS = {
    'jupyter': ['jupyter'],
    'matproj': ['atomate', 'pymongo'],
    'gpcam': ['gpcam'],
    'scrn-image': ['pygame', 'PIL'],
    'scrn-video': ['pygame', 'PIL', 'cv2'],
}

@functools.lru_cache(maxsize=128)
def _missing_dependency(extra):
    for module in EXTRAS[extra]:
        try:
            importlib.import_module(module)
        except ModuleNotFoundError:
            return module
    return False


@functools.lru_cache(maxsize=128)
def feature_loaded(extra):
    for module in EXTRAS[extra]:
        try:
            importlib.import_module(module)
        except ModuleNotFoundError:
            return False
    return True


def feature(extra, error=False):
    if extra not in EXTRAS:
        raise ValueError(f'{extra} is not a dtcs extra feature.')

    def requires_decorator(func):
        @functools.wraps(func)
        def decorated_func(*args, **kwargs):
            module = _missing_dependency(extra)
            if module:
                msg = f'Missing module \'{module}\'; '
                f'this function requires extra dependency \'{extra}\', '
                f'try `pip install dtcs[{extra}]`.'
                if error: raise ModuleNotFoundError(msg)
                else: warnings.warn(msg)
                return

            return func(*args, **kwargs)
        return decorated_func
    return requires_decorator


def flatten_dictionary(dic, prefix=tuple()):
    # Base case
    if not isinstance(dic, dict):
        return {prefix: dic}

    out_dict = {}
    for key, value in dic.items():
        sub_dict = flatten_dictionary(value, prefix=prefix + (key,))
        out_dict.update(sub_dict)
    return out_dict


def get_legend_patches(species):
    patch_info = {}
    for specie in species:
        color = color_map[specie]
        label = fr'${latex_map._get_value(specie, color=False)}$'
        patch_info[label] = color
    patches = [mpatches.Patch(color=color, label=label)
               for label, color in
               sorted(patch_info.items(), key=lambda x: x[0])]
    return patches
