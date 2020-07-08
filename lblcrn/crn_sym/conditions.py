"""Classes for describing the conditions of an reaction system experiment.

Exports:
    T: Time, accessible anywhere through sym.Symbol('t').
    ChemExpression: Abstract base class for classes in this module.
    Term: Apply a constant rate of change to a species.
    ConcEq: Set a species' concentration to be exactly this equation.
    ConcDiffEq: Set a species' rate of change to be exactly this equation.
    Conc: Define the initial concentration for a species.
    Schedule: Define amounts of a species to be added at a given time.

Usage:
    These classes are all created with the intention of passing it into a
    RxnSystem. For example:

    RxnSystem(
        ...
        Conc(x, 20)
        Schedule(y, {2: 30})
    )

    Use only one per species, as they don't logically collide (You can't fix
    two initial concentrations, that's a contradiction).
"""

from typing import Set

import monty.json
import sympy as sym
from sympy.parsing import sympy_parser

import lblcrn


T = sym.Symbol('t')  # Time


class ChemExpression(monty.json.MSONable):
    """Abstract base class for classes in this module.

    It is essentially a pair of a symbol/species and an expression which
    conveys information about that species.
    """

    def __init__(self, symbol: sym.Symbol, expression: sym.Expr):
        self.symbol = symbol
        self.expression = sym.sympify(expression)

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = {
            '@module': self.__class__.__module__,
            '@class': self.__class__.__name__,
            '@version': lblcrn.__version__,  # TODO: Better way to do this?
            'expression': str(self.expression),
            'symbol': str(self.symbol)
        }
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        d['expression'] = sympy_parser.parse_expr(d['expression'])
        d['symbol'] = sympy_parser.parse_expr(d['symbol'])
        return cls(**d)

    def get_symbols(self) -> Set[sym.Symbol]:
        """Get all of the symbols in the ChemExpression (except for time)."""

        symbols = self.symbol.free_symbols
        symbols.update(self.expression.free_symbols)
        symbols.discard(T)  # Time is not a species
        return symbols

    def __repr__(self):
        return f'{self.__class__.__name__}(symbol={repr(self.symbol)}, ' \
               f'expression={repr(self.expression)})'


class Term(ChemExpression):  # TODO(Andrew) Document here & beyond.
    """An additive term in the ODE for a symbol.

    If x' = ... + 2 and x' = ... + y*t are terms, then the ODE for x
    should look something like x' ... + 2 + y*t, for example.

    Incompatible with: ConcEq, ConcDiffEq
    """

    def __str__(self):
        return f'term: [{self.symbol}]\' = ... + {self.expression}'


class ConcEq(ChemExpression):
    """An equation for the concentration of a species.

    This guarantees that the concentration of the symbol will always be exactly
    equal to the expression, regardless of what any terms dictate.

    Incompatible with: Term, ConcDiffEq, Schedule, Conc
    """

    def __str__(self):
        return f'[{self.symbol}] = {self.expression}'


class ConcDiffEq(ChemExpression):
    """An equation for the derivative of the concentration of a species.

    This guarantees that the derivative of the concentration of the symbol will
    always be exactly equal to the expression, regardless of what any terms
    dictate.

    Incompatible with: CTerm, ConcEq,
    """

    def __str__(self):
        return f'[{self.symbol}]\' = {self.expression}'


class Schedule(monty.json.MSONable):
    """A schedule describing when amounts of a species are added/removed.

    Attributes:
        symbol: The species the Schedule corresponds to.
    """

    def __init__(self, symbol: sym.Symbol, schedule={}):  # TODO: Remove dict
        """Create a new schedule given a symbol and a schedule.

        The schedule can be either a dictionary or a list. Internally, it will
        keep the dictionary format.

        The dictionary format is {time: amount,}, where at time, it will add
        amount.

        The list format is [(time_difference, amount), ] where it will wait
        each time_difference and then add the amount. You can convert it to
        the dictionary format with a cumulative sum of time_differences.

        If you do not specify anything to do at time 0, it will assume that at
        time 0 you want concentration 0.
        """
        self.symbol = symbol

        # Handle schedule if it's a dict
        if isinstance(schedule, dict):
            self._schedule = schedule  # TODO: Get a better name for this

        # Handle schedule if it's a list or tuple
        elif isinstance(schedule, list):
            self._schedule = {}

            # Sum the times
            time = 0
            for time_diff, amount in schedule:
                time += time_diff
                self._schedule[time] = amount

        # Add the initial if it's not specified.
        if not 0 in self._schedule:
            self._schedule[0] = 0

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = {
            '@module': self.__class__.__module__,
            '@class': self.__class__.__name__,
            '@version': lblcrn.__version__,  # TODO: Better way to do this?
            'symbol': str(self.symbol),
            'schedule': self._schedule
        }
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        d['symbol'] = sympy_parser.parse_expr(d['symbol'])
        return cls(**d)

    def items(self):
        """So that you can do schedule.items() as if it's a dict."""
        return self._schedule.items()

    def __str__(self):
        s = f'{self.__class__.__name__} for {self.symbol}:\n'

        kvp = [pair for pair in self._schedule.items()]
        kvp.sort(key=lambda pair: pair[0])

        for time, amount in kvp:
            s += f'@t = {time} add {amount}\n'
        return s[:-1]

    def __repr__(self):
        return f'{self.__class__.__name__}(symbol={repr(self.symbol)}, ' \
               f'schedule={repr(self._schedule)})'


class Conc(Schedule):
    """Specify the initial concentration of a symbol.

    Internally, it is a schedule where the symbol gets its initial
    concentration added at t=0.
    """

    def __init__(self, symbol, concentration, locs=[]):
        """
        Make a new Conc: it's a Schedule where it all gets added at t=0.
        """

        Schedule.__init__(self, symbol, {0: concentration})
        self.concentration = concentration
        self.locs = locs

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict()
        del d['schedule']
        d['concentration'] = self.concentration
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        d['symbol'] = sympy_parser.parse_expr(d['symbol'])
        return cls(**d)

    def __str__(self):
        return f'[{self.symbol}] = {self.concentration} @ t=0'

    def __repr__(self):
        return f'{self.__class__.__name__}(symbol={repr(self.symbol)}, ' \
               f'concentration={repr(self.concentration)})'
