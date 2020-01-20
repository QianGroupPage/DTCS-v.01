"""
CRN Symulator:

This is a sympy-based chemistry interpretation engine.

It's essentially a really easy way to write chemical reactions in Python which are converted
into many different forms of initial value problem/ODE.

Note that it doesn't actually _solve_ the initial value problem. That is delegated to the solver.

Credits:
Dr. Jin Qian, Domas Buracas, Ye Wang, Andrew Bogdan, Rithvik Panchapakesan
"""

# *** Libraries ***
import numpy as np
import sympy as sym

from typing import List, Tuple, Set
import abc

from process_sympy_eqs import process_sympy_eqs

# *** Constants ***
T = sym.Symbol('t') # time

# *** Classes ***
class ChemExpression(abc.ABC):
    """
    A pair of a species and an expression (about that species)

    This class is an abstract superclass for structures of the form (species, expression), where
    the subclass gives meaning to the pair.
    """

    def __init__(self, species: sym.Symbol, expression: sym.Expr):
        self.species = species
        self.expression = sym.sympify(expression)

    def get_species(self) -> Set[sym.Symbol]:
        species = self.species.free_symbols
        species.update(self.expression.free_symbols)
        species.discard(T) # time is not a species
        return species

class Term(ChemExpression):
    """
    An additive term in the ODE for a species.

    If x' = ... + 2 and x' = ... + y*t are terms, then the ODE for x
    should look something like x' ... + 2 + y*t, for example.
    """

    def __str__(self):
        return 'term: [' + str(self.species) + ']\' = ... + ' + str(self.expression)
    def __repr__(self):
        return 'Term(species=' + repr(self.species) + ', expression=' + str(self.expression) + ')'

class ConcEq(ChemExpression):
    """
    An equation which describes the concentration of a species.

    This guarantees that the concentration of the species will always be exactly equal
    to the expression, regardless of what any terms dictate.
    """

    def __str__(self):
        return '[' + str(self.species) + '] = ' + str(self.expression)
    def __repr__(self):
        return 'ConcEq(species=' + repr(self.species) + ', expression=' + str(self.expression) + ')'

class ConcDiffEq(ChemExpression):
    """
    An equation which describes the derivative of the concentration of a species.

    This guarantees that the derivative of the concentration of the species will always 
    be exactly equal to the expression, regardless of what any terms dictate.
    """

    def __str__(self):
        return '[' + str(self.species) + ']\' = ' + str(self.expression)
    def __repr__(self):
        return 'ConcDiffEq(species=' + repr(self.species) + ', expression=' + str(self.expression) + ')'

class Rxn:
    """
    A chemical reaction with reactants, products, and a rate constant.
    """

    def __init__(self, reactants: sym.Expr, products: sym.Expr, k: float):
        """
        Create a new reaction by giving equation of the reactants, the equation of the products,
        and the rate constant.

        This is intended to look like reactants -> products @ rate k. That is, if
        your chemical equation looks like 2A -> B + C, then your reactants is 2A.
        That is, the reactants is not a list of reactants, but the equation on the left-hand
        side of the chemical equation. Similarly, products is the equation on the right-hand
        side of the chemical equation.
        """

        # Note that although the type suggestion is sym.Expr, that is a suggestion
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

    def get_species(self) -> Set[sym.Symbol]:
        species = set()
        species.update(self.reactants.free_symbols)
        species.update(self.products.free_symbols)
        return species

    def to_terms(self) -> List[Term]:
        """
        Create a list of terms from the reaction.

        Each term is essentially the reaction rate but positive or negative, 
        depending on if it is a reactand or the product.

        This uses the reaction rate formula r = k*(concentration of products^coefficient of product):
        if this is not true for your chemical equation, then you will get an error.
        """
        
        # Get the lefts and right 
        # Here we're assuming that reactants and products are linear, chemical equations.
        lefts = self.reactants.as_coefficients_dict()
        rights = self.products.as_coefficients_dict()

        # Make a dict (species : term), initilizing each to 0.
        term_dict = {}
        for species in self.get_species():
            term_dict[species] = sym.sympify(0)
        
        # Make the common part of all the terms
        common_part = self.rate_constant
        for species in self.reactants.free_symbols:
            common_part *= (species ** lefts[species])
        
        # Populate the lefts
        for species in self.reactants.free_symbols:
            term_dict[species] += -1 * lefts[species] * common_part
            
        # Populate the rights
        for species in self.products.free_symbols:
            term_dict[species] += rights[species] * common_part
            
        # Make it into terms
        terms = []
        for species in term_dict.keys():
            terms.append(Term(species, term_dict[species]))
            
        return terms

    def __str__(self):
        return str(self.reactants) + ' → ' + str(self.products) + ' @ k=' + str(self.rate_constant)
    def __repr__(self):
        return 'Rxn(reactants=' + repr(self.reactants) + ', products=' + repr(self.products) + ', k=' + str(self.rate_constant) + ')'

class RevRxn(Rxn):
    """
    A reversible reaction, essentially a reaction with two rate constants.

    Its use is to be quickly unpacked into two Rxns.
    """

    def __init__(self, reactants: sym.Expr, products: sym.Expr, k1: float, k2: float = None):
        """
        Create a reversible reaction by giving equation of the reactants, the equation of the products,
        and the rate constant.

        This is intended to look like reactants <-> products @ rate k1, with the reverse rate k2.
        If you don't specify a k2, it will assume that k2 = 1/k1.

        The reactants is not a list of reactants, but the equation on the left-hand
        side of the chemical equation. Similarly, products is the equation on the right-hand
        side of the chemical equation.
        """

        Rxn.__init__(self, reactants, products, k1)

        if k2 is None:
            self.rate_constant_reverse = 1 / k1
        else:
            self.rate_constant_reverse = k2

    def get_species(self) -> Set[sym.Symbol]:
        species = set()
        species.update(self.reactants.free_symbols)
        species.update(self.products.free_symbols)
        return species

    def to_rxns(self) -> Tuple[Rxn]:
        return Rxn(self.reactants, self.products, k=self.rate_constant), Rxn(self.products, self.reactants, k=self.rate_constant_reverse)

    def to_terms(self) -> List[Term]:
        rxns = self.to_rxns()
        return [*rxns[0].to_terms(), *rxns[1].to_terms()]
    
    def __str__(self):
        return str(self.reactants) + ' ↔ ' + str(self.products) + ' @ k1=' + str(self.rate_constant) + ', k2=' + str(self.rate_constant_reverse)

    def __repr__(self):
        return 'Revrxn(reactants=' + repr(self.reactants) + ', products=' + repr(self.products) + ', k1=' + str(self.rate_constant) + ', k2=' + str(self.rate_constant_reverse) + ')'

class Schedule:
    """
    A schedule describing when amounts of a species are added and removed.
    """

    def __init__(self, species: sym.Symbol, schedule={}):
        """
        Create a new schedule given a species and a description of the schedule.

        The scheduel can be either a dictionary or a list. Internally, it will keep
        the dictionary format.

        The dictionary format is {time: amount,}, where at time, it will add amount.

        The list format is [(time_differnece, amount), ] where it will wait each
        time_difference and then add the amount. You can convert it to the dictionary
        format with a cumulative sum of time_differences.

        If you do not specify anything to do at time 0, it will assume that at
        time 0 you want concentration 0.
        """

        self.species = species

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
        s = 'schedule: [' + str(self.species) + ']:\n'

        kvp = [pair for pair in self.schedule.items()]
        kvp.sort(key=lambda pair: pair[0])

        for time, amount in kvp:
            s += '@t = ' + str(time) + ' add ' + str(amount) + '\n'
        return s[:-1]

    def __repr__(self):
        return 'Schedule(species=' + repr(self.species) + ', schedule=' + repr(self.schedule) + ')'

class Conc(Schedule):
    """
    Specify the initial concentration of a species.

    Internally, it is a schedule where the species gets its initial concentration added at t=0.
    """

    def __init__(self, species, concentration):
        """
        Make a new Conc: it's a Schedule where it all gets added at t=0.
        """
        
        Schedule.__init__(self, species, {0: concentration})
        self.concentration = concentration

    def __str__(self):
        return '[' + str(self.species) + '] (@ t=0) = ' + str(self.concentration)
    def __repr__(self):
        return 'Conc(species=' + repr(self.species) + ', concentration=' + repr(self.concentration) + ')'

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
        self.components = flatter_components

        # Split into terms, schedules, and conc (diff.) eq.s
        self.terms = []
        self.schedules = []
        self.conc_eqs = []
        self.conc_diffeqs = []

        for component in self.components:
            if isinstance(component, Schedule):
                self.schedules.append(component)
            elif isinstance(component, Rxn):
                self.terms.extend(component.to_terms())
            elif isinstance(component, Term):
                self.terms.append(component)
            elif isinstance(component, ConcEq):
                self.conc_eqs.append(component)
            elif isinstance(component, ConcDiffEq):
                self.conc_diffeqs.append(component)
            else:
                assert False # !!! make an eeror?

        # Pick an order for the species
        self._species = set()
        for term in self.terms:
            self._species.update(term.get_species())
        for schedule in self.schedules:
            self._species.add(schedule.species)
        for equation in self.conc_eqs:
            self._species.update(equation.get_species())
        for equation in self.conc_diffeqs:
            self._species.update(equation.get_species())
        self._species = list(self._species)

        # Make an indexing dictionary
        self.species_index = {}
        for index, species in enumerate(self._species):
            self.species_index[species] = index

    def get_ode_expressions(self) -> List[sym.Expr]:
        """
        Return a list of expressions, corresponding to the derivative of the concentration of
        each species in the reaction system.

        The collection is ordered, and that order is accessible through species_index or
        through get_species.
        """
        
        # Make an emtpy ODE list
        odes = [sym.sympify(0)] * len(self._species)
        
        # Sum all terms for each species
        for term in self.terms:
            odes[self.species_index[term.species]] += term.expression

        # Set the conc diffeqs
        for equation in self.conc_diffeqs:
            odes[self.species_index[equation.species]] = equation.expression

        return odes

    def get_ode_functions(self):
        pass

    def get_species(self) -> List[sym.Symbol]:
        return self._species

    def __str__(self):
        s = 'rxn system with components:\n'
        for component in self.components:
            s += str(component) + '\n'
        return s[:-1]

    def __repr__(self):
        return 'RxnSystem(components=' + repr(self.components) + ')'

# *** Functions ***
def species(names: str) -> Tuple[sym.Symbol]:
    """A wrapper for sym.symbols"""
    return sym.symbols(names)
