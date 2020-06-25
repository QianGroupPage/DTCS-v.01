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


import sympy as sym
from typing import List, Set, Tuple

from lblcrn.crn_sym import conditions


class Rxn:  # TODO(Andrew) Document here & beyond.
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

    def get_symbols(self) -> Set[sym.Symbol]:
        symbol = set()
        symbol.update(self.reactants.free_symbols)
        symbol.update(self.products.free_symbols)
        return symbol

    # TODO: this is bulk crn-specific method, and as a result should be removed.
    def to_terms(self) -> List[conditions.Term]:
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

    def __str__(self):
        return str(self.reactants) + ' → ' + str(self.products) + ' @ k=' + str(self.rate_constant)

    def __repr__(self):
        return 'Rxn(reactants=' + repr(self.reactants) + ', products=' + repr(self.products) + ', k=' + \
               str(self.rate_constant) + ')'


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

    def get_symbols(self) -> Set[sym.Symbol]:
        symbol = set()
        symbol.update(self.reactants.free_symbols)
        symbol.update(self.products.free_symbols)
        return symbol

    def to_rxns(self) -> Tuple[Rxn, Rxn]:
        return Rxn(self.reactants, self.products, k=self.rate_constant), Rxn(self.products, self.reactants, k=self.rate_constant_reverse)

    def to_terms(self) -> List[conditions.Term]:
        rxns = self.to_rxns()
        return [*rxns[0].to_terms(), *rxns[1].to_terms()]
    
    def __str__(self):
        return str(self.reactants) + ' ↔ ' + str(self.products) + ' @ k1=' + str(self.rate_constant) + ', k2=' + str(self.rate_constant_reverse)

    def __repr__(self):
        return 'RevRxn(reactants=' + repr(self.reactants) + ', products=' + repr(self.products) + ', k1=' + str(self.rate_constant) + ', k2=' + str(self.rate_constant_reverse) + ')'
