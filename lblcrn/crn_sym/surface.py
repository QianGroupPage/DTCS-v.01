"""
CRN - surface.py:

Structures for the creation of surfaces

Credits:
Dr. Jin Qian, Domas Buracas, Andrew Bogdan, Rithvik Panchapakesan, Ye Wang
"""

# *** Libraries ***
from typing import List, Tuple, Union

# *** Classes ***
class Surface:
    """
    An orbital in a species, this is essentially a named tuple, it's a class for readability purposes.
    """

    def __init__(self, name: str, size: Tuple[int], color: Union[Tuple[int],List[int], str]=None):
        self.name = name
        self.size = size
        self.color = color

    def __str__(self):
        row = " ".join([self.name] * self.size[1])
        return "\n".join(row, self.size[0]) 

    def __repr__(self):
        return "Surface(name=" + self.name + ', size=' + str(self.size) + , ", color=" + str(self.color)')'
