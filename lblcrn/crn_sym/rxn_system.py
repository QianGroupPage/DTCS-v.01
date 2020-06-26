import copy
from typing import List
import sympy as sym

from lblcrn.crn_sym import species
from lblcrn.crn_sym import surface
from lblcrn.crn_sym import conditions
from lblcrn.crn_sym.reaction import Rxn
from lblcrn.crn_sym.surface_reaction import SurfaceRxn
from lblcrn.common import generate_new_color, color_to_RGB


class RxnSystem:
    """
    A collection of Terms, Rxns, Schedules, etc. which describes a chemical reaction system
    for the purpose of simulating it (namely the concentrations of each chemical) over time.
    """

    def __init__(self, *components):
        """
        Create a new reaction system by giving it Rxns, Revrxns, Concs, Schedules, Terms, ConcEqs, and
        ConcDiffEqs in any order.

        If you have a function returning a collection of the above, you do not have to worry about
        unpacking the collection: it will unpack and flatten lists and tuples for you.
        """

        # Flatten the components
        flat = False
        while not flat:
            flatter_components = []
            flat = True
            for component in components:
                if isinstance(component, list) or isinstance(component, tuple):
                    flat = False
                    flatter_components.extend(component)
                else:
                    flatter_components.append(component)
            components = flatter_components
        self.components = flatter_components  # pyint: disable=

        # Split into terms, schedules, and conc (diff.) eq.s
        # These are not sorted by the symbol index.
        self.terms = []
        self.schedules = []
        self.conc_eqs = []
        self.conc_diffeqs = []
        self.species_manager = None
        self.surface = None

        for component in self.components:
            if isinstance(component, conditions.Schedule):
                self.schedules.append(component)
            elif isinstance(component, SurfaceRxn):
                pass
            elif isinstance(component, Rxn):
                self.terms.extend(component.to_terms())
            elif isinstance(component, conditions.Term):
                self.terms.append(component)
            elif isinstance(component, conditions.ConcEq):
                self.conc_eqs.append(component)
            elif isinstance(component, conditions.ConcDiffEq):
                self.conc_diffeqs.append(component)
            elif isinstance(component, species.SpeciesManager):
                self.species_manager = component
            elif isinstance(component, surface.Surface):
                self.surface = component
            else:
                assert False, f'Unknown input {component} of type ' + str(type(component))

        # Pick an order for the symbol
        self._symbols = set()
        for term in self.terms:
            self._symbols.update(term.get_symbols())
        for schedule in self.schedules:
            self._symbols.add(schedule.symbol)
        for equation in self.conc_eqs:
            self._symbols.update(equation.get_symbols())
        for equation in self.conc_diffeqs:
            self._symbols.update(equation.get_symbols())
        self._symbols = list(self._symbols)

        # Make an indexing dictionary
        self.symbol_index = {}
        for index, symbol in enumerate(self._symbols):
            self.symbol_index[symbol] = index

        # Make symbol:concentration/scheudle list
        # Make default (Conc 0 @ t=0 for each species) scheduler list
        self.scheduler = [conditions.Conc(symbol, 0) for symbol in self._symbols]
        # Overwrite scheduler with Concs and Schedules
        for schedule in self.schedules:
            self.scheduler[self.symbol_index[schedule.symbol]] = schedule

        # an index from symbols to their colors
        self.color_index = None

    def get_ode_expressions(self) -> List[sym.Expr]:
        """
        Return a list of expressions, corresponding to the derivative of the concentration of
        each symbol in the reaction system.

        The collection is ordered, and that order is accessible through symbol_index or
        through get_symbols.
        """

        # Make an emtpy ODE list
        odes = [sym.sympify(0)] * len(self._symbols)

        # Sum all terms for each symbol
        for term in self.terms:
            odes[self.symbol_index[term.symbol]] += term.expression

        # Set the conc diffeqs
        for equation in self.conc_diffeqs:
            odes[self.symbol_index[equation.symbol]] = equation.expression

        return odes

    def get_ode_functions(self):
        """
        Return a function with signature func(t, y) representing the reaction system.
        """

        symbols = self.get_symbols()
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

    def get_conc_functions(self):
        """
        TODO
        Returns:
        """
        symbols = self.get_symbols()
        time = sym.symbols('t')

        conc_eq_funcs = {}
        for conc_eq in self.conc_eqs:
            index = self.symbol_index[conc_eq.symbol]
            func = sym.lambdify((time, symbols), conc_eq.expression)
            conc_eq_funcs[index] = func

        return conc_eq_funcs

    def get_species(self) -> List[species.Species]:
        species = []

        for symbol in self._symbols:
            species.append(self.species_manager[symbol])

        return species

    def get_symbols(self) -> List[sym.Symbol]:
        """
        Give all the symbols in the reaction system in a fixed order.
        """
        return copy.copy(self._symbols)

    # TODO
    def get_colors(self):
        """
        :return: colors for each species, if color is assigned.
        """
        colors = [] if not self.surface.color else [self.surface.color]
        if self.color_index is None:
            self.color_index = {}
        for index, symbol in enumerate(self._symbols):
            if self.color_index and symbol in self.color_index:
                continue
            color = self.species_manager[symbol].color
            if color is None:
                color = color_to_RGB(generate_new_color(colors))
                colors.append(color)
                self.species_manager[symbol].color = color

            self.color_index[symbol] = self.species_manager[symbol].color

        if self.surface:
            if self.color_index and self.surface.symbol in self.color_index:
                pass
            else:
                color = self.surface.color
                if color is None:
                    color = color_to_RGB(generate_new_color(colors))
                    colors.append(color)
                    self.surface.color = color
                self.color_index[self.surface.symbol()] = color
        return self.color_index




    def __str__(self):
        s = 'rxn system with components:\n'
        for component in self.components:
            s += str(component) + '\n'
        return s[:-1]

    def __repr__(self):
        return 'RxnSystem(components=' + repr(self.components) + ')'
