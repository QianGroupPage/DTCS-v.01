"""TODO(Andrew)"""

from __future__ import annotations

from typing import List

import copy

import monty.json
import sympy as sym

import lblcrn
from lblcrn.spec.xps import XPSSpecies as Species, XPSSpeciesManager as SpeciesManager
from lblcrn.spec.crn.bulk.rxn_system import BulkRxnSystem
from lblcrn.spec.crn.surface.reaction import SurfaceRxn, SurfaceRevRxn
from lblcrn.spec.crn.surface import surface
from lblcrn.spec.crn.bulk import conditions
from lblcrn.spec.model_input.state_variable import T, P

class SurfaceRxnSystem(BulkRxnSystem):

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
        for comp in components:
            if isinstance(comp, SpeciesManager):
                self.species_manager = comp
            elif isinstance(comp, surface.Surface):
                self.surface = comp
            else:
                clean_comps.append(comp)
        super().__init__(*clean_comps, name=name, description=description)
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
        """A list for names that appear on the surface.
        """
        return [self.surface.name] + [s.name for s in self.surface.sites]


    # def _update_elements(self, components):
    #     """
    #     Update the input elements in the Rsys.
    #     :return:
    #     """
    #     raise NotImplementedError()  # TODO(Andrew): Depreciate
    #
    #     assert len(components) > 0, 'Must pass at least one reaction and a species manager to a solution system.'
    #     if len(components) == 1:
    #         assert not isinstance(components[0], SpeciesManager), 'Must pass at least one reaction to a solution system.'
    #         raise AssertionError('Must pass a species manager to a solution system.')
    #
    #     # Flatten the components
    #     flat = False
    #     while not flat:
    #         flatter_components = []
    #         flat = True
    #         for component in components:
    #             if isinstance(component, list) or isinstance(component, tuple):
    #                 flat = False
    #                 flatter_components.extend(component)
    #             else:
    #                 flatter_components.append(component)
    #         components = flatter_components
    #     self.components = flatter_components  # pyint: disable=
    #
    #     # Split into terms, schedules, and conc (diff.) eq.s
    #     # These are not sorted by the symbol index.
    #     # self.terms = []  # Terms are terms.
    #     self.surface_rxns = []
    #     # self.schedules = []
    #     # self.conc_eqs = []
    #     # self.conc_diffeqs = []
    #     self.species_manager = None
    #     self.surface = None
    #     # self.reactions = []
    #     self.temperature = T(value=None)
    #     self.pressure = P(value=None)
    #
    #     for component in self.components:
    #     #     if isinstance(component, conditions.Schedule):
    #     #         self.schedules.append(component)
    #         if isinstance(component, Rxn):
    #             # self.terms.extend(component.to_terms())
    #             # pass
    #             self.surface_rxns.append(component)
    #     #     elif isinstance(component, Rxn):
    #     #         self.reactions.append(component)
    #     #         # TODO (ye): double check if this mismatches with Surface CRN capabilities
    #     #         self.terms.extend(component.to_terms())
    #     #     elif isinstance(component, conditions.Term):
    #     #         self.terms.append(component)
    #     #     elif isinstance(component, conditions.ConcEq):
    #     #         self.conc_eqs.append(component)
    #     #     elif isinstance(component, conditions.ConcDiffEq):
    #     #         self.conc_diffeqs.append(component)
    #         elif isinstance(component, surface.Surface):
    #             self.surface = component
    #         elif isinstance(component, SpeciesManager):
    #             self.species_manager = component
    #         elif isinstance(component, T):
    #             self.temperature = component
    #         elif isinstance(component, P):
    #             self.pressure = component
    #         # else:
    #         #     raise AssertionError(f'Unknown input {component} of type ' + str(type(component)) + ' to reaction system.')
    #
    #     assert self.species_manager is not None, 'Must pass a species manager to a solution system.'
    #
    #     # Share the surface name to the species manager
    #     if self.surface is not None:
    #         self.species_manager.default_surface_name = self.surface.name
    #
    #     self._update_symbols()
    #     self._generate_network_graph()
    #
    #     # Map the reactions to rates:
    #     self.rxns_to_rates = {rxn: [rxn.rate_constant, rxn.rate_constant_reverse] if isinstance(rxn, SurfaceRevRxn) else [rxn.rate_constant] for rxn in self.reactions}
    #     # self.rxns_by_name = {rxn.name: rxn for rxn in self.rxns}
    #     # TODO: support Surface CRN.
    #     self.surface_rxns_to_rates = {surface_rxn: surface_rxn.rate for surface_rxn in self.surface_rxns}
    #
    # def _update_symbols(self):
    #     """Iterate over reactions etc and update the symbols set and scheduler.
    #
    #     This method should be called whenever the reaction system is updated.
    #     """
    #     raise NotImplementedError()  # TODO(Andrew): Depreciate
    #
    #     # Pick an order for the symbol
    #     self._symbols = set()
    #     for rxn in self.surface_rxns:
    #         for s in rxn.get_symbols():
    #             if s not in self.surface.symbols:
    #                 self._symbols.add(s)
    #     for term in self.terms:
    #         self._symbols.update(term.get_symbols())
    #     for schedule in self.schedules:
    #         self._symbols.add(schedule.symbol)
    #     for equation in self.conc_eqs:
    #         self._symbols.update(equation.get_symbols())
    #     for equation in self.conc_diffeqs:
    #         self._symbols.update(equation.get_symbols())
    #     self._symbols = sorted(list(self._symbols), key=lambda s: str(s))
    #
    #     # Make an indexing dictionary
    #     self.symbol_index = {}
    #     for index, symbol in enumerate(self._symbols):
    #         self.symbol_index[symbol] = index
    #
    #     # Make symbol:concentration/schedule list
    #     # Make default (Conc 0 @ t=0 for each species) scheduler list
    #     self.scheduler = [conditions.Conc(symbol, 0) for symbol in self._symbols]
    #     # Overwrite scheduler with Concs and Schedules
    #     for schedule in self.schedules:
    #         self.scheduler[self.symbol_index[schedule.symbol]] = schedule
    #
    # # def get_ode_expressions(self) -> List[sym.Expr]:
    # #     """Return a list of expressions, corresponding to the derivative of the
    # #     concentration of each symbol in the reaction system.
    # #
    # #     The collection is ordered, and that order is accessible through
    # #     symbol_index or through get_symbols.
    # #     """
    # #
    # #     # Make an emtpy ODE list
    # #     odes = [sym.sympify(0)] * len(self._symbols)
    # #
    # #     # Sum all terms for each symbol
    # #     for term in self.terms:
    # #         odes[self.symbol_index[term.symbol]] += term.expression
    # #
    # #     # Set the conc diffeqs
    # #     for equation in self.conc_diffeqs:
    # #         odes[self.symbol_index[equation.symbol]] = equation.expression
    # #
    # #     return odes
    #
    # # def get_ode_functions(self):
    # #     """Return the ODE function with signature func(t, y) of the system."""
    # #     symbols = self.get_symbols()
    # #     odes = self.get_ode_expressions()
    # #     time = sym.symbols('t')
    # #     conc_eq_funcs = self.get_conc_functions()
    # #
    # #     # This makes a function with signature
    # #     # f((specie1, specie2, ...), time) -> (specie1, specie2, ...)
    # #     # Where 'specieN' is the concentration at that timestep.
    # #     # This is meant to be decorated and fed into SciPy's ODEINT solver.
    # #     undecorated_ode = sym.lambdify((time, symbols), odes)
    # #
    # #     # Decorate the function to add fixed conc equations
    # #     def decorated_ode(time: float, concs: List[float]):
    # #         for index, func in conc_eq_funcs.items():
    # #             concs[index] = func(time, concs)
    # #
    # #         return undecorated_ode(time, concs)
    # #
    # #     return decorated_ode
    #
    # # def get_conc_functions(self):
    # #     """Return the functions corresponding to the ConcEqs.
    # #
    # #     Used internally by get_ode_functions."""
    # #     symbols = self.get_symbols()
    # #     time = sym.symbols('t')
    # #
    # #     conc_eq_funcs = {}
    # #     for conc_eq in self.conc_eqs:
    # #         index = self.symbol_index[conc_eq.symbol]
    # #         func = sym.lambdify((time, symbols), conc_eq.expression)
    # #         conc_eq_funcs[index] = func
    # #     return conc_eq_funcs
    #
    # def _get_species(self) -> List[Species]:
    #     raise NotImplementedError()  # TODO(Andrew): Depreciate
    #     species = []
    #
    #     for symbol in self._symbols:
    #         species.append(self.species_manager[symbol])
    #
    #     return species
