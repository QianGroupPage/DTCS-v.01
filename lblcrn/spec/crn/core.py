"""TODO"""

from typing import Dict, List, Mapping, Set

import sympy as sym

from lblcrn.spec.crn.crn_abc import SymSpec
from lblcrn.spec.crn.rxn_system import RxnSystem
from lblcrn.spec.species import SpeciesCollection


class CRNSpec(SymSpec):

    _schema = [
        'rsys',
        'species',
        'sim_type',
    ]

    def __init__(self,
                 *components,
                 rsys: RxnSystem = None,
                 species: SpeciesCollection = None,
                 sim_type: str = None,
                 **kwargs):
        super().__init__(**kwargs)

        self.rsys: RxnSystem = rsys
        self.species: SpeciesCollection = species
        if sim_type:  # So as to not override any default
            self.sim_type = sim_type

        # If they give it components like it's a RxnSystem, deal with it
        rsys_components = []
        for component in components:
            if isinstance(component, SpeciesCollection):
                # If self.species, then they gave two SpeciesCollections.
                if self.species:
                    raise TypeError('You cannot give two SpeciesCollections.')
                self.species = component
            elif isinstance(component, RxnSystem):
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
            self.rsys = RxnSystem(*rsys_components)

    def get_symbols(self) -> Set[sym.Symbol]:
        return self.rsys.get_symbols()

    def get_symbols_ordered(self) -> List[sym.Symbol]:
        return self.rsys.get_symbols_ordered()

    @property
    def symbol_index(self) -> Dict[sym.Symbol, int]:
        return self.rsys.symbol_index

    def rename(self, mapping: Mapping):
        pass
        raise NotImplementedError()
