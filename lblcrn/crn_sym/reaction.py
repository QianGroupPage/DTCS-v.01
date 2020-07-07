"""Classes for defining a chemical reaction system.

Exports:
    Rxn: A chemical reaction
    RevRxn: A reversible chemical reaction
    RxnSystem: A collection of reactions and conditions (e.g. initial
        concentrations).

Usage:
    Usage essentially always looks like the following, where you make a
    reaction system with other classes as inputs.

    RxnSystem(
        sm,

        Rxn(x, y, k=3.2),
        RevRxn(x + y, 2z, k1=0.01, k2=100),

        Conc(x, 2),
        ...
    )

    It looks like this because it was made to imitate the predecessor
    mathematica project.
"""

import copy
from typing import List, Optional, Set, Tuple

import monty.json
import sympy as sym
from sympy.parsing import sympy_parser

import lblcrn
from lblcrn.crn_sym import species
from lblcrn.crn_sym import conditions


class Rxn(monty.json.MSONable):
    """A chemical reaction with reactants, products, and a rate constant.

    Attributes:
        reactants: A sym.Expr of the reactants.
        products: A sym.Expr of the products.
        rate_constant: A float, the rate constant of the chemical reaction.
    """

    def __init__(self, reactants: Optional[sym.Expr],
                 products: Optional[sym.Expr], k: float = 0):
        """Create a new reaction by giving equation of the reactants.

        This is intended to look like reactants -> products @ rate k. That is,
        if your chemical equation looks like 2A -> B + C, then your reactants
        is 2A.

        Args:
            reactants: The left-hand side of the chemical reaction.
            products: The right-hand side of the chemical reaction.
            k: The rate constant.
        """

        # Note that the type suggestion is Optional[sym.Expr].
        # It is possible that the user could pass None or 0 or 1 in.
        # Hence, sanitize input
        if not isinstance(reactants, sym.Expr):
            self.reactants = sym.sympify(0)
        else:
            self.reactants = reactants

        if not isinstance(products, sym.Expr):
            self.products = sym.sympify(0)
        else:
            self.products = products

        self.rate_constant = k

    def get_symbols(self) -> Set[sym.Symbol]:
        symbol = set()
        symbol.update(self.reactants.free_symbols)
        symbol.update(self.products.free_symbols)
        return symbol

    # TODO: this is bulk crn-specific method, and as a result should be removed.
    def to_terms(self) -> List[conditions.Term]:
        """Create a list of terms from the reaction.

        Each term is essentially the reaction rate but positive or negative, 
        depending on if it is a reactand or the product.

        This uses the reaction rate formula r =
        k*(concentration of products^coefficient of product):
        if this is false for your chemical equation, you will get an error.
        """
        
        # Get the lefts and right 
        # Here we're assuming that reactants and products are linear.
        lefts = self.reactants.as_coefficients_dict()
        rights = self.products.as_coefficients_dict()

        # Make a dict (symbol : term), initializing each to 0.
        term_dict = {}
        for symbol in self.get_symbols():
            term_dict[symbol] = sym.sympify(0)
        
        # Make the common part of all the terms
        common_part = self.rate_constant
        for symbol in self.reactants.free_symbols:
            common_part *= (symbol ** lefts[symbol])
        
        # Populate the lefts
        for symbol in self.reactants.free_symbols:
            term_dict[symbol] += -1 * lefts[symbol] * common_part
            
        # Populate the rights
        for symbol in self.products.free_symbols:
            term_dict[symbol] += rights[symbol] * common_part
            
        # Make it into terms
        terms = []
        for symbol in term_dict.keys():
            terms.append(conditions.Term(symbol, term_dict[symbol]))
            
        return terms

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = {
            '@module': self.__class__.__module__,
            '@class': self.__class__.__name__,
            '@version': lblcrn.__version__,  # TODO: Better way to do this?
            'reactants': str(self.reactants),
            'products': str(self.products),
            'k': self.rate_constant
        }
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        d['reactants'] = sympy_parser.parse_expr(d['reactants'])
        d['products'] = sympy_parser.parse_expr(d['products'])
        return cls(**d)


    def __str__(self):
        return f'{self.reactants} -> {self.products} @ k={self.rate_constant}'

    def __repr__(self):
        return f'{self.__class__.__name__}' \
               f'(reactants={repr(self.reactants)}, ' \
               f'products={repr(self.products)}, ' \
               f'k={self.rate_constant})'


class RevRxn(Rxn):
    """A reversible reaction, essentially a reaction with two rate constants.

    Its use is to be quickly unpacked into two Rxns.
    """

    def __init__(self, reactants: Optional[sym.Expr],
                 products: Optional[sym.Expr], k: float, k2: float = None):
        """Create a reversible reaction by giving equation.

        This is intended to look like reactants <-> products @ rate k1,
        with the reverse rate k2.

        Args:
            reactants: The left-hand side of the chemical reaction.
            products: The right-hand side of the chemical reaction.
            k: The rate constant.
            k2: Optional, the rate constant for the reverse reaction. If not
                supplied, it's assumed to be 1/k.
        """

        super().__init__(reactants, products, k)

        if k2 is None:
            self.rate_constant_reverse = 1 / k
        else:
            self.rate_constant_reverse = k2

    def get_symbols(self) -> Set[sym.Symbol]:
        symbol = set()
        symbol.update(self.reactants.free_symbols)
        symbol.update(self.products.free_symbols)
        return symbol

    def to_rxns(self) -> Tuple[Rxn, Rxn]:
        return Rxn(self.reactants, self.products, k=self.rate_constant), \
               Rxn(self.products, self.reactants, k=self.rate_constant_reverse)

    def to_terms(self) -> List[conditions.Term]:
        rxns = self.to_rxns()
        return [*rxns[0].to_terms(), *rxns[1].to_terms()]

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict()
        d['k2'] = self.rate_constant_reverse
        return d

    def __str__(self):
        return f'{self.reactants} <-> {self.products} ' \
               f'@ k={self.rate_constant}, k2={self.rate_constant_reverse}'

    def __repr__(self):
        return f'{self.__class__.__name__}' \
               f'(reactants={repr(self.reactants)}, ' \
               f'products={repr(self.products)}, ' \
               f'k={self.rate_constant}, k2={self.rate_constant_reverse})'


class RxnSystem(monty.json.MSONable):
    """A chemical reaction system, for simulation.

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

    def __init__(self, *components):
        """Create a new reaction system. Requires a SpeciesManager.

        Accepts Rxns, Revrxns, Concs, Schedules, Terms, ConcEqs,
        ConcDiffEqs, and (one) SpeciesManager in any order.

        If you have a function returning a collection of the above, you do
        not have to worry about unpacking the collection: it will unpack and
        flatten lists and tuples for you.
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
        self.components = flatter_components

        # Split into terms, schedules, and conc (diff.) eq.s
        # These are not sorted by the symbol index.
        self.terms = []
        self.schedules = []
        self.conc_eqs = []
        self.conc_diffeqs = []
        self.species_manager = None

        for component in self.components:
            if isinstance(component, conditions.Schedule):
                self.schedules.append(component)
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
            else:
                assert False, 'Unknown input type ' + str(type(component))

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

    def get_ode_expressions(self) -> List[sym.Expr]:
        """Return a list of expressions, corresponding to the derivative of the
        concentration of each symbol in the reaction system.

        The collection is ordered, and that order is accessible through
        symbol_index or through get_symbols.
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
        """Return the ODE function with signature func(t, y) of the system."""

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

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = {
            '@module': self.__class__.__module__,
            '@class': self.__class__.__name__,
            '@version': lblcrn.__version__,  # TODO: Better way to do this?
            'components': [comp.as_dict() for comp in self.components]
        }
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        decode = monty.json.MontyDecoder().process_decoded
        components = [decode(comp) for comp in d['components']]
        return cls(*components)

    def __str__(self):
        s = self.__class__.__name__ + ' with components:\n'
        for component in self.components:
            comp_lines = str(component).splitlines()
            s += ''.join([f'\t{line}\n' for line in comp_lines])
        return s[:-1]

    def __repr__(self):
        return f'{self.__class__.__name__}(components={repr(self.components)})'
