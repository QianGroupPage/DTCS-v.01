"""
CRN Symulator:

This is a sympy-based chemistry interpretation engine.

It's essentially a really easy way to write chemical reactions in Python which are converted
into many different forms of initial value problem/ODE.

Note that it doesn't actually _solve_ the initial value problem. That is delegated to the solver.

Credits:
Dr. Jin Qian, Domas Buracas, Ye Wang, Andrew Bogdan, Rithvik Panchapakesan
"""

from lblcrn.crn_sym.species import *
from lblcrn.crn_sym.conditions import *
from lblcrn.crn_sym.reaction import *
