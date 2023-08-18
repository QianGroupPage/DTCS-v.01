"""TODO"""

from typing import Dict, List, Mapping, Set
import copy
import warnings

import sympy as sym

from dtcs.spec.crn.sym_abc import SymSpec
from dtcs.spec.crn.rxn_abc import RxnSystemABC
from dtcs.spec.species import SpeciesManager
from dtcs.optim.gp import evaluate

from dtcs.optim.gibbs import CRNGibbsDataset
from dtcs.optim.sys_gen import system_generator, xps_generator


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
        self.sm: SpeciesManager = species

        # If they give it components like it's a RxnSystem, deal with it
        rsys_components = []
        for component in components:
            if isinstance(component, SpeciesManager):
                # If self.species, then they gave two SpeciesCollections.
                if self.sm:
                    raise TypeError('You cannot give two SpeciesCollections.')
                self.sm = component
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
                self.rsys.extend(rsys_components)
            else:
                self.rsys = self._rxn_sys_cls(*rsys_components)

    def get_symbols(self) -> Set[sym.Symbol]:
        return self.rsys.get_symbols()

    def get_symbols_ordered(self) -> List[sym.Symbol]:
        return self.rsys.get_symbols_ordered()

    def get_rates(self):
        return self.rsys.get_rates()

    def subs_gibbs(self, gibbs):
        """TODO(Andrew)"""
        crn = copy.deepcopy(self)
        crn.rsys = self.rsys.subs_gibbs(gibbs)
        return crn

    def subs_rates(self, rates):
        """TODO(Andrew)"""
        crn = copy.deepcopy(self)
        crn.rsys = self.rsys.subs_rates(rates)
        return crn

    def fit_rates(self, iterations, experimental, ignore):
        """TODO(Andrew)"""
        instrument = evaluate(
            crn=self,
            inst_args=dict(
                experimental=experimental,
                ignore=ignore,
            ),
            iterations=iterations,
        )
        return self.subs_rates(instrument.best_rates)

    def fit_gibbs(
            self,
            sample_at_ev,
            # gibbs_to_rates,
            num_energies,
            **system_kwargs
    ):
        sys = system_generator(
            crn=self,
            sample_at_ev=sample_at_ev,
            **system_kwargs,
        )

        dsg = CRNGibbsDataset.from_sim(
            sim=sys,
            crn=self,
            num_energies=num_energies,
        )

        # TODO(Andrew): This is a bodge for testing
        dsg._xps = xps_generator(
            crn=self,
            **system_kwargs,
        )

        return dsg

    @property
    def symbol_index(self) -> Dict[sym.Symbol, int]:
        return self.rsys.symbol_index

    def rename(self, mapping: Mapping):
        pass
        raise NotImplementedError()
