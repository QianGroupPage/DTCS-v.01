"""TODO(Andrew)"""

import copy

import sympy as sym

from lblcrn.spec.crn.sym_abc import ChemInfo


class Coverage(ChemInfo):
    """TODO(Andrew)"""

    def __init__(self, symbol: sym.Symbol, coverage: float, **kwargs):
        super().__init__(symbol=symbol, **kwargs)

        self.coverage = coverage
        if not (0 <= coverage <= 1):
            raise ValueError('Coverage should be between 0 and 1.')
