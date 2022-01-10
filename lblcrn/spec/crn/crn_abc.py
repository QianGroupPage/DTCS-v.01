"""TODO"""

from typing import Dict, List, Mapping, Set
import copy

import sympy as sym

from lblcrn.spec.crn.sym_abc import SymSpec
from lblcrn.spec.crn.rxn_abc import RxnSystemABC
from lblcrn.spec.species import SpeciesManager


# TODO(Andrew) move this to abc
class CRNSpecABC(SymSpec):

    _rxn_sys_cls = RxnSystemABC

    def __init__(self,
                 *components,
                 rsys: RxnSystemABC = None,
                 species: SpeciesManager = None,
                 sim_type: str = None,  # TODO(Andrew) Remove
                 **kwargs):
        super().__init__(**kwargs)

        self.rsys: RxnSystemABC = rsys
        self.species: SpeciesManager = species
        if sim_type:  # So as to not override any default
            self.sim_type = sim_type

        # If they give it components like it's a RxnSystem, deal with it
        rsys_components = []
        for component in components:
            if isinstance(component, SpeciesManager):
                # If self.species, then they gave two SpeciesCollections.
                if self.species:
                    raise TypeError('You cannot give two SpeciesCollections.')
                self.species = component
            elif isinstance(component, self._rxn_sys_cls):
                # If self.rsys, then they gave two RxnSystems.
                if self.rsys:
                    raise TypeError('You cannot give two RxnSystems.')
                self.rsys = component
            else:
                # We're going to create a reaction system with these
                rsys_components.append(component)

        if rsys_components:
            if self.rsys:
                # They gave a rsys, so we can't create another
                raise TypeError('Attempting to create a RxnSystem given these '
                                'inputs, but you already supplied one.')
            self.rsys = self._rxn_sys_cls(*rsys_components)

    def get_symbols(self) -> Set[sym.Symbol]:
        return self.rsys.get_symbols()

    def get_symbols_ordered(self) -> List[sym.Symbol]:
        return self.rsys.get_symbols_ordered()

    def get_rates(self):
        return self.rsys.get_rates()

    def subs_rates(self, rates):
        """TODO(Andrew)"""
        crn = copy.deepcopy(self)
        crn.rsys = self.rsys.subs_rates(rates)
        return crn

    def calc_rates(self):
        raise NotImplementedError()

    @property
    def symbol_index(self) -> Dict[sym.Symbol, int]:
        return self.rsys.symbol_index

    def rename(self, mapping: Mapping):
        pass
        raise NotImplementedError()
