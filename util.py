from typing import Dict
from sympy import lambdify, parse_expr

SympySymbol = str

def sympy_to_list(odes: Dict[SympySymbol, str]):

