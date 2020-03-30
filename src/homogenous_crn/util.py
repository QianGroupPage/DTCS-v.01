from typing import Dict
from sympy import lambdify
import re

SympySymbol = str

# def sympy_to_list(odes: Dict[SympySymbol, str]):


def multiple_replace(dict, text):
  """
  :param dict: a dictionary from the string to be replaced to what it should be replaced to
  :return: text with correct replacements.
  """
  # Create a regular expression  from the dictionary keys
  regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))
  # For each match, look-up corresponding value in dictionary
  return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text)

