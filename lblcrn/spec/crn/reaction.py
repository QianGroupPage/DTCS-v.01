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

from typing import List, Mapping, Optional, Set, Tuple

from monty.json import jsanitize
import sympy as sym
from sympy.parsing import sympy_parser

from lblcrn.spec.crn.crn_abc import SymSpec


class Rxn(SymSpec):
    """A chemical reaction with reactants, products, and a rate constant.

    Attributes:
        reactants: A sym.Expr of the reactants.
        products: A sym.Expr of the products.
        rate_constant: A float, the rate constant of the chemical reaction.
    """

    _schema = [
        'reactants',
        'products',
        'rate_constant',
    ]

    def __init__(self, reactants: Optional[sym.Expr],
                 products: Optional[sym.Expr], k: float, **kwargs):
        """Create a new reaction by giving equation of the reactants.

        This is intended to look like reactants -> products @ rate k. That is,
        if your chemical equation looks like 2A -> B + C, then your reactants
        is 2A.

        Args:
            reactants: The left-hand side of the chemical reaction.
            products: The right-hand side of the chemical reaction.
            k: The rate constant.
        """
        super().__init__(**kwargs)
        self.rate_constant = k

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

    def get_symbols(self) -> Set[sym.Symbol]:
        symbols = set()
        symbols.update(self.reactants.free_symbols)
        symbols.update(self.products.free_symbols)
        return symbols

    def rename(self, mapping: Mapping):
        self.reactants.subs(mapping)
        self.products.subs(mapping)

    def as_dict(self, sanitize=True) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict(sanitize=False)
        d['reactants'] = str(d['reactants'])
        d['products'] = str(d['products'])
        if sanitize:
            d = jsanitize(d, strict=True)
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        d['reactants'] = sympy_parser.parse_expr(d['reactants'])
        d['products'] = sympy_parser.parse_expr(d['products'])
        d['k'] = d.pop('rate_constant')
        return super(Rxn, cls).from_dict(d)

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

    _schema = [
        'rate_constant_reverse',
    ]

    def __init__(self, reactants: Optional[sym.Expr],
                 products: Optional[sym.Expr], k: float,
                 k2: Optional[float] = None, **kwargs):
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

        super().__init__(reactants=reactants, products=products, k=k, **kwargs)
        self.rate_constant_reverse = k2 or 1 / self.rate_constant

    def to_rxns(self) -> Tuple[Rxn, Rxn]:
        return Rxn(self.reactants, self.products, k=self.rate_constant), \
               Rxn(self.products, self.reactants, k=self.rate_constant_reverse)

    @classmethod
    def from_dict(cls, d: dict):
        d['k2'] = d.pop('rate_constant_reverse')
        return super(RevRxn, cls).from_dict(d)

    def __str__(self):
        return f'{self.reactants} <-> {self.products} ' \
               f'@ k={self.rate_constant}, k2={self.rate_constant_reverse}'

    def __repr__(self):
        return f'{self.__class__.__name__}' \
               f'(reactants={repr(self.reactants)}, ' \
               f'products={repr(self.products)}, ' \
               f'k={self.rate_constant}, k2={self.rate_constant_reverse})'
