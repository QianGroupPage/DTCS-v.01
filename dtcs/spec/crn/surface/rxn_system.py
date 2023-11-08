"""TODO(Andrew)"""

from __future__ import annotations

from typing import List

import copy

import monty.json
import sympy as sym

from dtcs.spec.xps import XPSSpecies as Species, XPSSpeciesManager as SpeciesManager
from dtcs.spec.crn.rxn_abc import RxnSystemABC
from dtcs.spec.crn.surface.reaction import SurfaceRxn, SurfaceRevRxn
from dtcs.spec.crn.surface import surface
from dtcs.spec.crn.surface.conditions import Coverage
from dtcs.spec.crn.bulk.rxn_system import BulkRxnSystem

class SurfaceRxnSystem(RxnSystemABC):

    """A chemical reaction system, for simulation.

    A collection of Terms, Rxns, Schedules, etc. which describes a chemical
    reaction system for the purpose of simulating it (namely the concentrations
    of each chemical) over time.

    :var components: Everything the RxnSystem contains
    :var terms: Terms in the ODE of the system.
    :var reactions: Bulk CRN reactions in the system.
    :var schedules: The Schedules and Concs passed during initialization.
    :var conc_eqs: The ConcEqs in the system.
    :var conc_diffeqs: The ConcDiffEqs in the system.
    :var species_manager: The SpeciesManager the system uses.
    :var symbol_index: A dictionary {sym.Symbol: int} to keep track of the order of the symbols.
    :var scheduler: A comprehensive Schedule of the system, has entries (which might be Conc(species, 0) for each species which isn't set by a ConcEq.
    :var surface_names: A list for names that appear on the surface.
    :var network_graph: A graph representing the reaction network
    :var network_graph_pos: A position map for the network graph.
    """

    def __init__(self, *components, name: str = "", description: str = ""):
        """Create a new reaction system. Requires a SpeciesManager.

        Accepts Rxns, Revrxns, Concs, Schedules, Terms, ConcEqs,
        ConcDiffEqs, and (one) SpeciesManager in any order.

        If you have a function returning a collection of the above, you do
        not have to worry about unpacking the collection: it will unpack and
        flatten lists and tuples for you.
        """
        clean_comps = []
        elements = []
        for comp in components:
            if isinstance(comp, SpeciesManager):
                self.species_manager = comp
            elif isinstance(comp, surface.Surface):
                self.surface = comp
            elif isinstance(comp, Coverage):
                elements.append(comp)
            else:
                clean_comps.append(comp)
        super().__init__(*clean_comps, name=name)
        self.elements.extend(elements)
        # A partial list of all properties, with initial value.
        # Relation to govern (T, P) and the rate for each species.
        self.tprate_relation = None
        # Relation to govern (T, P) and the gas concentration.
        self.tpconc_relation = None
        #self._update_elements(components)

        #self._update_symbols()
        #self._generate_network_graph()

        # Map the reactions to rates:
        self.rxns_to_rates = {rxn: [rxn.rate_constant, rxn.rate_constant_reverse] if isinstance(rxn, SurfaceRevRxn) else [rxn.rate_constant] for rxn in self.reactions}
        # self.rxns_by_name = {rxn.name: rxn for rxn in self.rxns}
        # TODO: support Surface CRN.
        self.surface_rxns_to_rates = {surface_rxn: surface_rxn.rate for surface_rxn in self.surface_rxns}

    @property
    def surface_rxns(self) -> List[SurfaceRxn]:
        return self.by_subclass()[SurfaceRxn]

    @property
    def components(self):  # TODO(Andrew) I am very rushed, sorry
        return self.elements

    @property
    def surface_names(self) -> List[str]:
        """A list for names that appear on the surface."""
        return [self.surface.name] + [s.name for s in self.surface.sites]

    def to_bulk(self):
        """
        TODO(Andrew)
        """
        reactions = []
        for elem in self.elements:
            if isinstance(elem, SurfaceRxn):
                reactions.append(elem.to_bulk())

        return BulkRxnSystem(*reactions)