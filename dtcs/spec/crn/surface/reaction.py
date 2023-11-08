
from typing import Mapping, Optional, Sequence, Set, Tuple, Union
import sympy as sym

from dtcs.spec.crn.rxn_abc import RxnABC, RevRxnABC, Relation
from dtcs.spec.crn.bulk.reaction import BulkRxn, BulkRevRxn
from dtcs.common import util, display


class SurfaceRxn(RxnABC):
    """A chemical reaction with reactants, products, and a rate constant,
    tailored to map with the physical interactions between species.
    """

    def __init__(
            self,
            reactants: Union[Tuple[sym.Symbol], sym.Symbol],
            products: Union[Tuple[sym.Symbol], sym.Symbol],
            k: Optional[Relation] = None,
            dg: Optional[Relation] = None,
            **kwargs
    ):
        """Create a new reaction by giving a tuple of reactants, a tuple of
        the products, and the rate constant.

        This intends to represent that:
        (Species 1 at location A, Species 2 at adjacent location B) ->
        (Species 3 at location A, Species 4 at location B).
        """
        super().__init__(
            reactants=None,
            products=None,
            k=k,
            dg=dg,
            **kwargs,
        )

        # --- Class Variables ---
        self.reactants: Tuple
        self.products: Tuple

        # --- Digest Reactants & Products ---
        # The type hint is Optional[Union[Tuple[sym.Symbol], sym.Symbol]].
        # Make everything into a Sequence type, as order matters
        if not isinstance(reactants, Sequence): reactants = (reactants, )
        if not isinstance(products, Sequence): products = (products, )

        # --- Error Checking ---
        # This could be done with list comprehension, but the longer code makes
        #  for more verbose error output.
        self.reactants = []
        for reactant in reactants:
            if not isinstance(reactant, sym.Symbol):
                raise TypeError(f'Must supply Symbols, not {type(reactant)}')
            self.reactants.append(reactant)

        self.products = []
        for product in products:
            if not isinstance(product, sym.Symbol):
                raise TypeError(f'Must supply Symbols, not {type(product)}')
            self.products.append(product)

        if len(self.reactants) != len(self.products):
            raise ValueError(f'Reactants and products must be equal length.')

        self.reactants = tuple(self.reactants)
        self.products = tuple(self.products)

    # --- Abstract Method Implementations -------------------------------------
    # TODO: reactants are input as symbols and need sm to get the species object.
    def get_symbols(self) -> Set[sym.Symbol]:
        symbols = super().get_symbols()
        symbols.update(self.reactants)
        symbols.update(self.products)
        return symbols

    def rename(self, mapping: Mapping):
        super().rename(mapping)
        self.reactants = (react.subs(mapping) for react in self.reactants)
        self.products = (prod.subs(mapping) for prod in self.products)

    def _latex_rxn(self) -> str:
        react_str = ', '.join((display.pretty_sym_subs(react) for react in self.reactants))
        prod_str = ', '.join((display.pretty_sym_subs(prod) for prod in self.products))
        return f'[{react_str}]' + r' \longrightarrow ' + f'[{prod_str}]'

    # --- Other Methods -------------------------------------------------------
    def to_bulk(self):
        """TODO(Andrew)"""
        return BulkRxn(sum(self.reactants),
                       sum(self.products),
                       k=self.rate_constant)

    # TODO: produce a version that takes in the species manager, replace
    #  default sites with top sites on appropriate occasions

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

    @util.depreciate
    def is_adsorption(self, sm, gas=None):
        """
        Rule: a reaction is adsorption if all reactants are gas or surface.
        :return:
        """
        if gas is None:
            return all([sm.is_gas(r) or sm.is_surface(r) for r in self.reactant_symbols])
        return all([(sm.is_gas(r) and r.name == gas) or sm.is_surface(r) for r in self.reactant_symbols])

    @util.depreciate
    def is_desorption(self, sm, gas=None):
        """
        Rule: a reaction is desorption if all products are gas or surface.
        """
        if gas is None:
            return all([sm.is_gas(r) or sm.is_surface(r) for r in self.product_symbols])
        else:
            return all([(sm.is_gas(r) and r.name == gas) or sm.is_surface(r) for r in self.product_symbols])

    @util.depreciate
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

    @property
    @util.depreciate
    def is_reversible(self):
        return False

    # def __str__(self):
    #     return self.reactants_str + ' → ' + self.products_str + ' @ k=' + str(self.rate_constant)
    #
    # def __repr__(self):
    #     return 'SurfaceRxn(reactants=' + repr(self.reactants) + ', products=' + repr(self.products) + ', k=' + \
    #            str(self.rate_constant) + ')'


class SurfaceRevRxn(SurfaceRxn, RevRxnABC):
    """
    A reversible surface reaction, essentially a surface reaction with two rate constants.

    Its use is to be quickly unpacked into two SurfaceRxns.
    """

    _rxn_cls = SurfaceRxn

    # def __init__(
    #     self,
    #     reactants: Union[Tuple[sym.Symbol], sym.Symbol],
    #     products: Union[Tuple[sym.Symbol], sym.Symbol],
    #     k: Optional[Relation] = None,
    #     k2: Optional[Relation] = None,
    #     dg: Optional[Relation] = None,
    #     **kwargs
    # ):
    #     """
    #     Create a reversible surface reaction by giving tuple of the reactants, the tuple of the products,
    #     and the rate constant.
    #
    #     This is intended to look like reactants <-> products @ rate k1, with the reverse rate k2.
    #     If you don't specify a k2, it will assume that k2 = 1/k1.
    #     """
    #
    #     SurfaceRxn.__init__(self, reactants, products, k=k)
    #
    #     if k2 is None:
    #         self.rate_constant_reverse = 1 / k1
    #     else:
    #         self.rate_constant_reverse = k2
    #
    #     self.forward_surface_rxn = SurfaceRxn(self.reactants, self.products, self.rate_constant)
    #     self.backward_surface_rxn = SurfaceRxn(self.products, self.reactants, k=self.rate_constant_reverse)
    #     self.is_reversible = True

    # --- Abstract Method Implementations -------------------------------------
    def _latex_rxn(self) -> str:
        react_str = ', '.join((display.pretty_sym_subs(react) for react in self.reactants))
        prod_str = ', '.join((display.pretty_sym_subs(prod) for prod in self.products))
        return f'[{react_str}]' + r' \longleftrightarrow ' + f'[{prod_str}]'

    # --- Other Methods -------------------------------------------------------
    def to_bulk(self):
        """TODO(Andrew)"""
        return BulkRevRxn(sum(self.reactants),
                          sum(self.products),
                          k=self.rate_constant,
                          k2=self.rate_constant_reverse)

    # def to_rxns(self) -> Tuple[SurfaceRxn, SurfaceRxn]:
    #     return self.forward_surface_rxn.surface_engine_str, self.backward_surface_rxn.surface_engine_str

    @property
    def surface_engine_str(self):
        """
        Mainly used to match between this class, and Surface CRN engine's internal rule class.
        :return: a str representation in exact same style output by the Surface CRN engine.
        """
        return self.reactants_str + " --> " + self.products_str

    @property
    @util.depreciate
    def is_reversible(self):
        return True

    # def __str__(self):
    #     return self.reactants_str + ' ↔ ' + self.products_str + ' @ k1=' + str(self.rate_constant) + ', ' \
    #            + 'k2=' + str(self.rate_constant_reverse)
    #
    # def __repr__(self):
    #     return 'SurfaceRevRxn(reactants=' + repr(self.reactants) + ', products=' + repr(self.products) + ', k1=' + str(
    #         self.rate_constant) + ', k2=' + str(self.rate_constant_reverse) + ')'
