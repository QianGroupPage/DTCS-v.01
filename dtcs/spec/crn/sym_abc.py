import abc
import copy
from typing import Set, Mapping

from monty.json import jsanitize
import sympy as sym
from sympy.parsing import sympy_parser

from dtcs.common.const import T
from dtcs.spec.spec_abc import Spec


class SymSpec(Spec):
    """TODO"""

    @abc.abstractmethod
    def get_symbols(self) -> Set[sym.Symbol]:
        """Return all the sym.Symbols in the spec."""
        pass

    @abc.abstractmethod
    def rename(self, mapping: Mapping):
        """Replace sym.Symbols with other sym.Symbols."""
        pass

    #@abc.abstractmethod
    def name(self):
        """SymSpecs can get generate their names from their data"""
        return str(self)

    def subs(self, mapping: Mapping):
        """Return a copy with substituted symbols."""
        clone = copy.copy(self)
        clone.rename(mapping)
        return clone


class ChemInfo(SymSpec):
    """TODO"""

    def __init__(self, symbol: sym.Symbol, **kwargs):
        super().__init__(**kwargs)
        self.symbol: sym.Symbol = symbol

    def get_symbols(self) -> Set[sym.Symbol]:
        return {self.symbol}

    @property
    def species(self):
        return str(self.symbol)

    def rename(self, mapping: Mapping):
        self.symbol.subs(mapping)

    def as_dict(self, sanitize=True) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict(sanitize=False)
        d['symbol'] = str(d['symbol'])
        if sanitize:
            d = jsanitize(d)
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        d['symbol'] = sympy_parser.parse_expr(d['symbol'])
        return super(ChemInfo, cls).from_dict(d)


class ChemExpression(ChemInfo):
    """Abstract base class for classes in this module.

    It is essentially a pair of a symbol/species and an expression which
    conveys information about that species.

    TODO: Say time-dependent
    """

    def __init__(self, symbol: sym.Symbol, expression: sym.Expr, **kwargs):
        super().__init__(symbol=symbol, **kwargs)
        self.expression: sym.Expr = sym.sympify(expression)

    def get_symbols(self) -> Set[sym.Symbol]:
        """Get all of the symbols in the ChemExpression (except for time)."""
        symbols = super().get_symbols()
        symbols.update(self.expression.free_symbols)
        symbols.discard(T)  # Time is not a species
        return symbols

    def rename(self, mapping: Mapping):
        """Replace sym.Symbols with other sym.Symbols."""
        super().rename(mapping)
        self.expression = self.expression.subs(mapping)

    def as_dict(self, sanitize=True) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict(sanitize=False)
        d['expression'] = str(d['expression'])
        if sanitize:
            d = jsanitize(d, strict=True)
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        d['expression'] = sympy_parser.parse_expr(d['expression'])
        return super(ChemExpression, cls).from_dict(d)

    def __repr__(self):
        return f'{self.__class__.__name__}(symbol={repr(self.symbol)}, ' \
               f'expression={repr(self.expression)})'
