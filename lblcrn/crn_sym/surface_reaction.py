from lblcrn.crn_sym.reaction import Rxn
from lblcrn.crn_sym.surface import Surface, Site
from typing import Set, Tuple
import sympy as sym


class SurfaceRxn(Rxn):
    """
    A chemical reaction with reactants, products, and a rate constant, tailored for to map with the physical'
    interactions between species.
    """

    def __init__(self, reactants: Tuple[sym.Symbol], products: Tuple[sym.Symbol], k: float):
        """
        Create a new reaction by giving a tuple of reactants, a tuple of the products,
        and the rate constant.

        This intends to represent that: (Species 1 at location A, Species 2 at adjacent location B)
        -> (Species 3 at location A, Species 4 at location B).
        """
        if len(reactants) != len(products):
            raise Exception(f"left={reactants}, right={products}, left and right hands of a surface reaction must"
                            + " be equal length.")

        for s in reactants + products:
            if s is not None and (not isinstance(s, sym.Symbol) and not isinstance(s, Surface)) and not isinstance(s, Site):
                raise Exception(f"{s} is not a of sympy.Symbol, Site, or Surface class. \n" +
                                "please create it using the sp method of a species manager")

        self.reactants = []
        for s in reactants:
            # if s is not None and (not isinstance(s, sym.Symbol) and not isinstance(s, Surface)):
            #     raise Exception(f"{s} is not a of sympy.Symbol class or Surface class. \n" +
            #                     "please create it using the sp method of a species manager")
            if isinstance(s, Surface):
                self.reactants.append(s.symbol())
            elif isinstance(s, Site):
                self.reactants.append(s.symbol)
            else:
                self.reactants.append(s)

        self.products = []
        for s in products:
            if isinstance(s, Surface):
                self.products.append(s.symbol())
            elif isinstance(s, Site):
                self.products.append(s.symbol)
            else:
                self.products.append(s)

        self.reactants = tuple(self.reactants)
        self.products = tuple(self.products)
        self.rate_constant = k

    # TODO: reactants are input as symbols and need sm to get the species object.
    def get_symbols(self) -> Set[sym.Symbol]:
        symbol = set()
        symbol.update(self.reactants)
        symbol.update(self.products)
        return symbol

    def to_terms(self) -> None:
        raise Exception("SurfaceRxn doesn't support to_terms method.")

    @property
    def reactants_str(self):
        return " + ".join([s.name for s in self.reactants])

    @property
    def products_str(self):
        return " + ".join([s.name for s in self.products])

    def __str__(self):
        return self.reactants_str + ' → ' + self.products_str + ' @ k=' + str(self.rate_constant)

    def __repr__(self):
        return 'SurfaceRxn(reactants=' + repr(self.reactants) + ', products=' + repr(self.products) + ', k=' + \
               str(self.rate_constant) + ')'


class SurfaceRevRxn(SurfaceRxn):
    """
    A reversible surface reaction, essentially a surface reaction with two rate constants.

    Its use is to be quickly unpacked into two SurfaceRxns.
    """

    def __init__(self, reactants: Tuple[sym.Symbol], products: Tuple[sym.Symbol], k1: float, k2: float = None):
        """
        Create a reversible surface reaction by giving tuple of the reactants, the tuple of the products,
        and the rate constant.

        This is intended to look like reactants <-> products @ rate k1, with the reverse rate k2.
        If you don't specify a k2, it will assume that k2 = 1/k1.
        """

        SurfaceRxn.__init__(self, reactants, products, k1)

        if k2 is None:
            self.rate_constant_reverse = 1 / k1
        else:
            self.rate_constant_reverse = k2

    def to_rxns(self) -> Tuple[Rxn, Rxn]:
        return SurfaceRxn(self.reactants, self.products, self.rate_constant), \
               SurfaceRxn(self.products, self.reactants, k=self.rate_constant_reverse)

    def __str__(self):
        return self.reactants_str + ' ↔ ' + self.products_str + ' @ k1=' + str(self.rate_constant) + ', ' \
               + 'k2=' + str(self.rate_constant_reverse)

    def __repr__(self):
        return 'SurfaceRevRxn(reactants=' + repr(self.reactants) + ', products=' + repr(self.products) + ', k1=' + str(
            self.rate_constant) + ', k2=' + str(self.rate_constant_reverse) + ')'
