import warnings
from typing import List, Callable, Dict

import sympy as sym

from lblcrn.spec.crn.bulk.conditions import ConcEq, ConcDiffEq, Term, Schedule, Conc
from lblcrn.spec.crn.bulk.reaction import BulkRxn
from lblcrn.spec.crn.rxn_system import RxnSystem
from lblcrn.spec.crn.core import CRNSpec
from lblcrn.spec.crn.reaction import Rxn


class BulkRxnSystem(RxnSystem):
    """TODO: Old!

    A chemical reaction system, for simulation.

    A collection of Terms, Rxns, Schedules, etc. which describes a chemical
    reaction system for the purpose of simulating it (namely the concentrations
    of each chemical) over time.

    Attributes:
        components: Everything the RxnSystem contains
        terms: Terms in the ODE of the system.
        schedules: The Schedules and Concs passed during initialization.
        conc_eqs: The ConcEqs in the system.
        conc_diffeqs: The ConcDiffEqs in the system.
        species_manager: The SpeciesManager the system uses.
        symbol_index: A dictionary {sym.Symbol: int} to keep track of the order
            of the symbols.
        scheduler: A comprehensive Schedule of the system, has entries (which
            might be Conc(species, 0) for each species which isn't set by a
            ConcEq.
    """

    @property
    def conc_eqs(self) -> List[ConcEq]:
        return list(filter(lambda c: isinstance(c, ConcEq),
                           self.components))

    @property
    def conc_diffeqs(self) -> List[ConcDiffEq]:
        return list(filter(lambda c: isinstance(c, ConcDiffEq),
                           self.components))

    @property
    def terms(self) -> List[Term]:
        terms = []
        for component in self.components:
            if isinstance(component, Rxn):
                if not isinstance(component, BulkRxn):
                    warnings.warn(f'{component} is a Rxn but not a BulkRxn, '
                                  f'please use BulkRxn when specifying bulk'
                                  f'CRN simulations.')
                # TODO: Perhaps do a try-catch and throw an explicit
                #  SpecError or something like it, if this errs?
                terms.extend(BulkRxn.to_terms(component))
        return terms

    @property
    def schedules(self) -> List[Schedule]:
        return list(filter(lambda c: isinstance(c, Schedule),
                           self.components))

    @property
    def scheduler(self) -> List[Schedule]:
        scheduler = [Conc(symbol, 0) for symbol in self.get_symbols_ordered()]
        for schedule in self.schedules:
            scheduler[self.symbol_index[schedule.symbol]] = schedule
        return scheduler

    def get_ode_expressions(self) -> List[sym.Expr]:
        """Return a list of expressions, corresponding to the derivative of the
        concentration of each symbol in the reaction system.

        The collection is ordered, and that order is accessible through
        symbol_index or through get_symbols.
        """

        # Make an emtpy ODE list
        odes = [sym.sympify(0)] * len(self.get_symbols())

        # Sum all terms for each symbol
        for term in self.terms:
            odes[self.symbol_index[term.symbol]] += term.expression

        # Set the conc diffeqs
        for equation in self.conc_diffeqs:
            odes[self.symbol_index[equation.symbol]] = equation.expression

        return odes

    def get_ode_functions(self) -> Callable:
        """Return the ODE function with signature func(t, y) of the system."""

        symbols = self.get_symbols_ordered()
        odes = self.get_ode_expressions()
        time = sym.symbols('t')
        conc_eq_funcs = self.get_conc_functions()

        # This makes a function with signature
        # f((specie1, specie2, ...), time) -> (specie1, specie2, ...)
        # Where 'specieN' is the concentration at that timestep.
        # This is meant to be decorated and fed into SciPy's ODEINT solver.
        undecorated_ode = sym.lambdify((time, symbols), odes)

        # Decorate the function to add fixed conc equations
        def decorated_ode(time: float, concs: List[float]):
            for index, func in conc_eq_funcs.items():
                concs[index] = func(time, concs)

            return undecorated_ode(time, concs)

        return decorated_ode

    def get_conc_functions(self) -> Dict[int, Callable]:
        """Return the functions corresponding to the ConcEqs.

        Used internally by get_ode_functions."""
        symbols = self.get_symbols()
        time = sym.symbols('t')

        conc_eq_funcs = {}
        for conc_eq in self.conc_eqs:
            index = self.symbol_index[conc_eq.symbol]
            func = sym.lambdify((time, symbols), conc_eq.expression)
            conc_eq_funcs[index] = func

        return conc_eq_funcs
