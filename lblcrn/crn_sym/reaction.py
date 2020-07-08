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
