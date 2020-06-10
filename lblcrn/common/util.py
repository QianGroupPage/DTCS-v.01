import re

import numpy as np


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