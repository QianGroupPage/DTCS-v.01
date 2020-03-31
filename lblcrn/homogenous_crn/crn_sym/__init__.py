"""
CRN Symulator:

This is a sympy-based chemistry interpretation engine.

It's essentially a really easy way to write chemical reactions in Python which are converted
into many different forms of initial value problem/ODE.

Note that it doesn't actually _solve_ the initial value problem. That is delegated to the solver.

Credits:
Dr. Jin Qian, Domas Buracas, Ye Wang, Andrew Bogdan, Rithvik Panchapakesan
"""

from .species import *
from .conditions import *
from .reaction import *
from .bindingenergy import *