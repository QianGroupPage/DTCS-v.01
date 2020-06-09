"""
...
"""

# *** Libraries ***
import sympy as sym
import numpy as np

import abc
import copy

from typing import List, Tuple, Union, Set

# *** Constants ***
T = sym.Symbol('t') # time

# *** Classes ***
class ChemExpression(abc.ABC):
    """
    A pair of a symbol and an expression (about that symbol)

    This class is an abstract superclass for structures of the form (symbol, expression), where
    the subclass gives meaning to the pair.
    """

    def __init__(self, symbol: sym.Symbol, expression: sym.Expr):
        self.symbol = symbol
        self.expression = sym.sympify(expression)

    def get_symbols(self) -> Set[sym.Symbol]:
        symbol = self.symbol.free_symbols
        symbol.update(self.expression.free_symbols)
        symbol.discard(T) # time is not a symbol
        return symbol


class Term(ChemExpression):
    """
    An additive term in the ODE for a symbol.

    If x' = ... + 2 and x' = ... + y*t are terms, then the ODE for x
    should look something like x' ... + 2 + y*t, for example.
    """

    def __str__(self):
        return 'term: [' + str(self.symbol) + ']\' = ... + ' + str(self.expression)
    def __repr__(self):
        return 'Term(symbol=' + repr(self.symbol) + ', expression=' + str(self.expression) + ')'


class ConcEq(ChemExpression):
    """
    An equation which describes the concentration of a symbol.

    This guarantees that the concentration of the symbol will always be exactly equal
    to the expression, regardless of what any terms dictate.
    """

    def __str__(self):
        return '[' + str(self.symbol) + '] = ' + str(self.expression)
    def __repr__(self):
        return 'ConcEq(symbol=' + repr(self.symbol) + ', expression=' + str(self.expression) + ')'


class ConcDiffEq(ChemExpression):
    """
    An equation which describes the derivative of the concentration of a symbol.

    # TODO: are these really initial value?
    This guarantees that the derivative of the concentration of the symbol will always 
    be exactly equal to the expression, regardless of what any terms dictate.
    """

    def __str__(self):
        return '[' + str(self.symbol) + ']\' = ' + str(self.expression)
    def __repr__(self):
        return 'ConcDiffEq(symbol=' + repr(self.symbol) + ', expression=' + str(self.expression) + ')'


class Schedule:
    """
    A schedule describing when amounts of a symbol are added and removed.
    """

    def __init__(self, symbol: sym.Symbol, schedule={}):
        """
        Create a new schedule given a symbol and a description of the schedule.

        The schedule can be either a dictionary or a list. Internally, it will keep
        the dictionary format.

        The dictionary format is {time: amount,}, where at time, it will add amount.

        The list format is [(time_difference, amount), ] where it will wait each
        time_difference and then add the amount. You can convert it to the dictionary
        format with a cumulative sum of time_differences.

        If you do not specify anything to do at time 0, it will assume that at
        time 0 you want concentration 0.
        """

        self.symbol = symbol

        # Handle schedule if it's a dict
        if isinstance(schedule, dict):
            self.schedule = schedule

        # Handle schedule if it's a list or tuple
        elif isinstance(schedule, list):
            self.schedule = {}

            # Sum the times
            time = 0
            for time_diff, amount in schedule:
                time += time_diff
                self.schedule[time] = amount

        # Add the initial if it's not specified.
        if not 0 in self.schedule:
            self.schedule[0] = 0

    def __str__(self):
        s = 'schedule: [' + str(self.symbol) + ']:\n'

        kvp = [pair for pair in self.schedule.items()]
        kvp.sort(key=lambda pair: pair[0])

        for time, amount in kvp:
            s += '@t = ' + str(time) + ' add ' + str(amount) + '\n'
        return s[:-1]

    def __repr__(self):
        return 'Schedule(symbol=' + repr(self.symbol) + ', schedule=' + repr(self.schedule) + ')'


class Conc(Schedule):
    """
    Specify the initial concentration of a symbol.

    Internally, it is a schedule where the symbol gets its initial concentration added at t=0.
    """

    def __init__(self, symbol, concentration):
        """
        Make a new Conc: it's a Schedule where it all gets added at t=0.
        """
        
        Schedule.__init__(self, symbol, {0: concentration})
        self.concentration = concentration

    def __str__(self):
        return '[' + str(self.symbol) + '] (@ t=0) = ' + str(self.concentration)
    def __repr__(self):
        return 'Conc(symbol=' + repr(self.symbol) + ', concentration=' + repr(self.concentration) + ')'