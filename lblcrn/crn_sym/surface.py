"""
CRN - surface.py:

Structures for the creation of surfaces

Credits:
Dr. Jin Qian, Domas Buracas, Andrew Bogdan, Rithvik Panchapakesan, Ye Wang
"""

# *** Libraries ***
from typing import List, Tuple, Union
import sympy as sym
from lblcrn.common import color_to_RGB

# *** Classes ***
class Surface:
    """
    A surface structure, by default, it is square with only top sites.
    """

    def __init__(self, name: str, size: Tuple[int], color: Union[Tuple[int], List[int], str] = None):
        self.name = name
        self.size = size

        if color:
            self.color = color_to_RGB(color)
        else:
            self.color = color

    #property
    def symbol(self):
        return sym.Symbol(self.name)

    def __str__(self):
        row = " ".join([self.name] * self.size[1])
        if self.color:
            return "\n".join([row] * self.size[0]) + f"\n color={self.color}\n" 
        return "\n".join([row] * self.size[0]) 

    def __repr__(self):
        return "Surface(name=" + self.name + ', size=' + repr(self.size) + " , color=" + repr(self.color) + ')'
