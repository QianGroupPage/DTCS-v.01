from dtcs.spec.crn.bulk.reaction import BulkRxn, BulkRevRxn
from dtcs.spec.crn.surface.surface import Surface, Site
from dtcs.spec.crn.surface.marker import Marker
from typing import Set, Tuple
import sympy as sym


class Markers():
    """
    A wrapper for the wrappers associated with a SurfaceRxn
    """
    def __init__(self):
        self.reactants = []
        self.products = []


class SurfaceRxn(BulkRxn):
    """
    A chemical reaction with reactants, products, and a rate constant, tailored to map with the physical
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
            if s is not None and (not isinstance(s, sym.Symbol) and not isinstance(s, Surface)) \
                    and not isinstance(s, Site) and not isinstance(s, Marker):
                raise Exception(f"{s} is not an object of sympy.Symbol, Site, Surface, or Marker class. \n" +
                                "please create it using the sp method of the species manager")
        self.markers = Markers()
        self.markers.reactants = []
        self.markers.products = []

        self.reactants = []
        for s in reactants:
            # if s is not None and (not isinstance(s, sym.Symbol) and not isinstance(s, Surface)):
            #     raise Exception(f"{s} is not a of sympy.Symbol class or Surface class. \n" +
            #                     "please create it using the sp method of a species manager")
            if isinstance(s, Surface):
                # TODO(Andrew): Maybe assume that it's the first site, if you give the whole surface
                self.reactants.append(s.symbol())
            elif isinstance(s, Site):
                self.reactants.append(s.symbol)
            # Handle a marker
            elif isinstance(s, Marker):
                self.reactants.append(s.species_symbol)
                self.markers.reactants.append(s)
            else:
                self.reactants.append(s)

        self.products = []
        for s in products:
            if isinstance(s, Surface):
                self.products.append(s.symbol())
            elif isinstance(s, Site):
                self.products.append(s.symbol)
            elif isinstance(s, Marker):
                self.products.append(s.species_symbol)
                self.markers.products.append(s)
            else:
                self.products.append(s)

        self.reactants = tuple(self.reactants)
        self.products = tuple(self.products)
        self.rate_constant = k                                                  #property rate_constant of surfacerevrxn object has no setter
        self.is_reversible = False

    # TODO: reactants are input as symbols and need sm to get the species object.
    def get_symbols(self) -> Set[sym.Symbol]:
        symbol = set()
        symbol.update(self.reactants)
        symbol.update(self.products)
        return symbol

    def to_bulk(self):
        """TODO(Andrew)"""
        return BulkRxn(sum(self.reactants),
                       sum(self.products),
                       k=self.rate_constant)

    # TODO: produce a version that takes in the species manager, replace default sites with top sites on appropriate
    # occasions

    @property
    def reactants_str(self):
        return " + ".join([s.name for s in self.reactants])

    @property
    def products_str(self):
        return " + ".join([s.name for s in self.products])

    @property
    def surface_engine_str(self):
        """
        Mainly used to match between this class, and Surface CRN engine's internal rule class.
        :return: a str representation in exact same style output by the Surface CRN engine.
        """
        return self.reactants_str + " --> " + self.products_str

    @rate_constant.setter
    def rate_constant(self, k):
        self.rate_constant = k

    def is_adsorption(self, sm, gas=None):
        """
        Rule: a reaction is adsorption if all reactants are gas or surface.
        :return:
        """
        if gas is None:
            return all([sm.is_gas(r) or sm.is_surface(r) for r in self.reactant_symbols])
        return all([(sm.is_gas(r) and r.name == gas) or sm.is_surface(r) for r in self.reactant_symbols])

    def is_desorption(self, sm, gas=None):
        """
        Rule: a reaction is desorption if all products are gas or surface.
        """
        if gas is None:
            return all([sm.is_gas(r) or sm.is_surface(r) for r in self.product_symbols])
        else:
            return all([(sm.is_gas(r) and r.name == gas) or sm.is_surface(r) for r in self.product_symbols])

    def sorption_species(self, sm):
        """
        Find the adsorption and desorption species.

        :param sm:
        :return:
        """
        if self.is_adsorption(sm=sm):
            return [r for r in self.reactant_symbols if sm.is_gas(r)][0]
        if self.is_desorption(sm=sm):
            return [r for r in self.product_symbols if sm.is_gas(r)][0]

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

        self.forward_surface_rxn = SurfaceRxn(self.reactants, self.products, self.rate_constant)
        self.backward_surface_rxn = SurfaceRxn(self.products, self.reactants, k=self.rate_constant_reverse)
        self.is_reversible = True

    def to_bulk(self):
        """TODO(Andrew)"""
        return BulkRevRxn(sum(self.reactants),
                          sum(self.products),
                          k=self.rate_constant,
                          k2=self.rate_constant_reverse)

    def to_rxns(self) -> Tuple[SurfaceRxn, SurfaceRxn]:
        return self.forward_surface_rxn.surface_engine_str, self.backward_surface_rxn.surface_engine_str

    @property
    def surface_engine_str(self):
        """
        Mainly used to match between this class, and Surface CRN engine's internal rule class.
        :return: a str representation in exact same style output by the Surface CRN engine.
        """
        return self.reactants_str + " --> " + self.products_str

    def __str__(self):
        return self.reactants_str + ' ↔ ' + self.products_str + ' @ k1=' + str(self.rate_constant) + ', ' \
               + 'k2=' + str(self.rate_constant_reverse)

    def __repr__(self):
        return 'SurfaceRevRxn(reactants=' + repr(self.reactants) + ', products=' + repr(self.products) + ', k1=' + str(
            self.rate_constant) + ', k2=' + str(self.rate_constant_reverse) + ')'
