"""TODO: Old

Classes for describing the conditions of an reaction system experiment.

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

import copy

import sympy as sym

from lblcrn.spec.crn.sym_abc import ChemInfo, ChemExpression


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


class Schedule(ChemInfo):
    """A schedule describing when amounts of a species are added/removed.

    Attributes:
        symbol: The species the Schedule corresponds to.
    """

    _default = {
        'schedule': dict,
    }

    def __init__(self, symbol: sym.Symbol, schedule=None, **kwargs):
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
        super().__init__(symbol=symbol, **kwargs)

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

        # # Add the initial if it's not specified.
        # if 0 not in self.schedule:
        #     self.schedule[0] = 0

    def unpack(self):
        return self.schedule.items()

    @classmethod
    def from_dict(cls, d: dict):
        if 'schedule' in d:
            d['schedule'] = {float(time): amount for time, amount
                             in d['schedule'].items()}
        return super(Schedule, cls).from_dict(d)

    def __str__(self):
        s = f'{self.__class__.__name__} for {self.symbol}:\n'

        kvp = [pair for pair in self.schedule.items()]
        kvp.sort(key=lambda pair: pair[0])

        for time, amount in kvp:
            s += f'@t = {time} add {amount}\n'
        return s[:-1]

    def __repr__(self):
        return f'{self.__class__.__name__}(symbol={repr(self.symbol)}, ' \
               f'schedule={repr(self.schedule)})'

    # TODO(Andrew) Copied from old:
    def update_conc(self, func=lambda x: x, inplace=False):
        """
        Update each concentration value in the schedule with a function of that value

        :param func: function to update concentration value; by default, it's concentration.
        :param inplace: if true, modify this object; otherwise, modify other objects;
        :return: None if inplace, otherwise a deepcopy of current schedule with modified values.
        """
        if inplace:
            for k in self.schedule:
                self.schedule[k] = func(self.schedule[k])
        else:
            new_schedule = copy.deepcopy(self)
            new_schedule = func(new_schedule)
            return new_schedule

    def items(self):
        """So that you can do schedule.items() as if it's a dict."""
        return self.schedule.items()

    @property
    def initial_concentration(self) -> float:
        """
        Find the initial concentration;
        :return: the initial concentration.
        """
        return self.schedule[0]


class Conc(Schedule):
    """Specify the initial concentration of a symbol.

    Internally, it is a schedule where the symbol gets its initial
    concentration added at t=0.
    """

    def __init__(self, symbol: sym.Symbol, concentration: float, **kwargs):
        """Make a new Conc: it's a Schedule where it all gets added at t=0."""
        super().__init__(symbol=symbol, schedule={0: concentration}, **kwargs)

    @property
    def concentration(self) -> float:
        return self.schedule[0]

    @concentration.setter
    def concentration(self, value: float):
        self.schedule[0] = value

    @concentration.deleter
    def concentration(self):
        del self.schedule

    # TODO: Change this to modify spec, so that spec has concentration in it

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        d['concentration'] = d.pop('schedule')['0']
        return super(Conc, cls).from_dict(d)

    def __str__(self):
        return f'[{self.symbol}] = {self.concentration} @ t=0'

    def __repr__(self):
        return f'{self.__class__.__name__}(symbol={repr(self.symbol)}, ' \
               f'concentration={repr(self.concentration)})'