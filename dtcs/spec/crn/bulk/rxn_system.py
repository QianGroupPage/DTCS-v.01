import copy
import random
import warnings
from collections import defaultdict
from typing import Tuple, List, Dict, Optional, Callable

import monty.json
# import networkx as nx
import sympy as sym
import numpy as np

#from jupyter_dash import JupyterDash
#import dash_cytoscape as cyto
#import dash_html_components as html
#import dash_core_components as dcc
#import dash

import dtcs
from dtcs.common import util
from dtcs.common.colors.color_gradient import color_to_HEX, color_to_RGB
from dtcs.spec.crn.bulk import conditions
from dtcs.spec.crn.bulk.conditions import ConcEq, ConcDiffEq, Term, Schedule, Conc
from dtcs.spec.crn.bulk.reaction import BulkRxn, BulkRevRxn
from dtcs.spec.crn.rxn_abc import RxnABC, RxnSystemABC
from dtcs.spec.spec_abc import SpecCollection
from dtcs.spec.model_input.state_variable import T, P
from dtcs.spec.model_input.relations.tprate_relation import TPRateRelation
from dtcs.spec.model_input.relations.ideal_gas_law import IdealGasLawRelation


# TODO(Andrew) Make this RxnSystem as an alias and RxnSystem -> RxnSystemABC
class BulkRxnSystem(RxnSystemABC):
    """TODO: Old!

    A chemical reaction system, for simulation.

    A collection of Terms, Rxns, Schedules, etc. which describes a chemical
    reaction system for the purpose of simulating it (namely the concentrations
    of each chemical) over time.

    Attributes:
        components: Everything the RxnSystem contains
        terms: Terms in the ODE of the system.
        schedules: The Schedules and Concs passed during initialization.
        conc_eqs: The ConcEqs in the system.
        conc_diffeqs: The ConcDiffEqs in the system.
        species_manager: The SpeciesManager the system uses.
        symbol_index: A dictionary {sym.Symbol: int} to keep track of the order
            of the symbols.
        scheduler: A comprehensive Schedule of the system, has entries (which
            might be Conc(species, 0) for each species which isn't set by a
            ConcEq.
    """
    @property
    def conc_eqs(self) -> List[ConcEq]:
        return self.by_subclass()[ConcEq]

    @property
    def conc_diffeqs(self) -> List[ConcDiffEq]:
        return self.by_subclass()[ConcDiffEq]

    @property
    def terms(self) -> List[Term]:
        terms = self.by_subclass()[Term]
        terms.extend([rxn.to_terms() for rxn in self.by_subclass()[BulkRxn]])

        # for component in self.components:
        #     if isinstance(component, RxnABC):
        #         if not isinstance(component, BulkRxn):
        #             warnings.warn(f'{component} is a Rxn but not a BulkRxn, '
        #                           f'please use BulkRxn when specifying bulk'
        #                           f'CRN simulations.')
        #         # TODO: Perhaps do a try-catch and throw an explicit
        #         #  SpecError or something like it, if this errs?
        #         terms.extend(BulkRxn.to_terms(component))
        return list(util.flat(terms))

    @property
    def schedules(self) -> List[Schedule]:
        return self.by_subclass()[Schedule]

    # @property
    # def color_index(self):
    #     return {index: self.elements[index].color for index in range(len(self))}

    @property
    @util.depreciate
    def scheduler(self) -> List[Schedule]:
        scheduler = [Conc(symbol, 0) for symbol in self.species_symbols]
        for schedule in self.schedules:
            scheduler[self.symbol_index[schedule.symbol]] = schedule
        return scheduler

    @property
    def schedule(self) -> Dict[float, Dict[str, float]]:
        scheduler = defaultdict(lambda: defaultdict(float))
        scheduler[0] = defaultdict(float)

        for schedule in self.schedules:
            for time, amount in schedule.schedule.items():
                scheduler[time][schedule.species] = amount

        return dict(scheduler)

    def get_ode_expressions(self) -> List[sym.Expr]:
        """Return a list of expressions, corresponding to the derivative of the
        concentration of each symbol in the reaction system.

        The collection is ordered, and that order is accessible through
        symbol_index or through get_symbols.
        """

        # Make an emtpy ODE list
        # TODO: Get symbols ordered???
        odes = [sym.sympify(0)] * len(self.get_symbols())

        # Sum all terms for each symbol
        for term in self.terms:
            odes[self.symbol_index[term.symbol]] += term.expression

        # Set the conc diffeqs
        for equation in self.conc_diffeqs:
            odes[self.symbol_index[equation.symbol]] = equation.expression

        return odes

    def get_ode_functions(self) -> Callable:
        """Return the ODE function with signature func(t, y) of the system."""

        symbols = self.get_symbols_ordered()
        odes = self.get_ode_expressions()
        time = sym.symbols('t')
        conc_eq_funcs = self.get_conc_functions()

        # This makes a function with signature
        # f((specie1, specie2, ...), time) -> (specie1, specie2, ...)
        # Where 'specieN' is the concentration at that timestep.
        # This is meant to be decorated and fed into SciPy's ODEINT solver.
        undecorated_ode = sym.lambdify((time, symbols), odes)

        # Decorate the function to add fixed conc equations
        def decorated_ode(time: float, concs: List[float]):
            for index, func in conc_eq_funcs.items():
                concs[index] = func(time, concs)

            return undecorated_ode(time, concs)

        return decorated_ode

    def get_conc_functions(self) -> Dict[int, Callable]:
        """Return the functions corresponding to the ConcEqs.

        Used internally by get_ode_functions."""
        symbols = self.get_symbols()
        time = sym.symbols('t')

        conc_eq_funcs = {}
        for conc_eq in self.conc_eqs:
            index = self.symbol_index[conc_eq.symbol]
            func = sym.lambdify((time, symbols), conc_eq.expression)
            conc_eq_funcs[index] = func

        return conc_eq_funcs

    # TODO(Andrew) look at
    @util.depreciate
    def get_colors(self):
        """Assign (if applicable) and return colors for all species.
        """
        raise NotImplementedError()

        # return {sym: next(_COLORS) for sym in self.get_symbols()}

        if self.color_index:
            return self.color_index

        random.seed(3)
        colors = []
        if self.surface and self.surface.color:
            color = color_to_RGB(self.surface.color)
            colors.append(color)

        if self.color_index is None:
            self.color_index = {}
        for index, symbol in enumerate(self.species_manager.symbols_ordering):
            if self.surface and symbol in self.surface.symbols:
                continue
            if self.color_index and symbol in self.color_index:
                continue
            color = self.species_manager[symbol].color
            if color is None:
                color = color_to_RGB(generate_new_color(colors))
                colors.append(color)
                self.species_manager[symbol].color = color
            else:
                color = color_to_RGB(color)

            self.color_index[symbol] = color

        if self.surface:
            if self.surface.color is None:
                color = color_to_RGB(generate_new_color(colors))
                colors.append(color)
                self.surface.color = color
            self.color_index[self.surface.symbol()] = self.surface.color

            for s in self.surface.sites:
                if not s.color:
                    color = color_to_RGB(generate_new_color(colors))
                    s.color = color
                else:
                    color = color_to_RGB(s.color)
                colors.append(color)
                self.color_index[s.symbol] = color

            for marker_name in self.species_manager.get_marker_names():
                color = color_to_RGB(generate_new_color(colors))
                marker_colors = set()
                for marker in self.species_manager.get_markers(marker_name):
                    if not marker.color:
                        marker.color = color
                    else:
                        color = color_to_RGB(marker.color)
                        marker_colors.add(color)
                if len(marker_colors) > 1:
                    raise ValueError(f"Marker with name {marker_name} was assigned multiple colors: " +
                                     ", ".join(marker_colors))

                colors.append(color)
                # TODO: using a string as key here, whereas all other keys are symbol.symbols
                self.color_index[marker_name] = color
        return self.color_index

    # TODO(Andrew) Uncomment when I install the dependencies
    # def _add_to_network_graph(self, r, p, w):
    #     """Add the given reactions and products to the network graph.
    #     """
    #     for reactant in r:
    #         for product in p:
    #             self.network_graph.add_edge(str(product), str(reactant), weight=w)
    #
    # def _generate_network_graph(self):
    #     """Create a reaction network graph (data structure).
    #
    #     Each species becomes a vertex and reactions are represented with directed edges from
    #     reactants to products (where the reaction constants are the edge weights).
    #     """
    #     self.network_graph = nx.DiGraph()
    #
    #     for rxn in self.reactions:
    #         r = rxn.reactants.free_symbols
    #         p = rxn.products.free_symbols
    #
    #         self._add_to_network_graph(r, p, rxn.rate_constant)
    #         if isinstance(rxn, RevRxn):
    #             self._add_to_network_graph(p, r, rxn.rate_constant_reverse)
    #     self.network_graph_pos = nx.spring_layout(self.network_graph)
    #     for k, pos in self.network_graph_pos.items():
    #         # Scale positions to make graphing easier.
    #         self.network_graph_pos[k] = np.array([500*pos[0], 500*pos[1]])
    #
    # def plot(self):
    #     """Plot the reaction network as an interactive graph.
    #
    #     The reaction network graph is plotted in a user-draggable view. Users can also add new
    #     reactions to the system using the provided text input (simply re-run the cell to view the
    #     updated graph with the new reaction).
    #     """
    #     G = self.network_graph
    #     pos = self.network_graph_pos
    #     largest_weight = -1
    #
    #     # Convert the networkx graph to a format that Cytoscape can use
    #     for node in G.nodes():
    #         for e in G.out_edges(node):
    #             weight = G.get_edge_data(node, e[1])["weight"]
    #             largest_weight = max(largest_weight, weight)
    #
    #     elements = []
    #     for node in G.nodes():
    #         x, y = pos[node]
    #         elements.append({"data": {"id": str(node), "label": str(node)}, "position": {"x": x, "y": y}})
    #         for e in G.out_edges(node):
    #             weight = G.get_edge_data(node, e[1])["weight"]
    #             elements.append({"data": {
    #                 "source": str(node), "target": str(e[1]), "weight": weight, "normalized_weight": weight / largest_weight * 3.75 + 2
    #             }})
    #
    #     # Display the plot as well as the related edit inputs
    #     cy = cyto.Cytoscape(
    #         id="network-plot",
    #         layout={"name": "preset"},
    #         style={"width": "100%", "height": "500px"},
    #         elements=elements,
    #         stylesheet=[
    #             {
    #                 "selector": "edge",
    #                 "style": {
    #                     # The default curve style does not work with certain arrows
    #                     "curve-style": "bezier",
    #                     "source-arrow-color": "black",
    #                     "source-arrow-shape": "triangle",
    #                     "label": "data(weight)",
    #                     "width": "data(normalized_weight)",
    #                 }
    #             },
    #             {
    #                 "selector": "node",
    #                 "style": {
    #                     "label": "data(label)",
    #                 }
    #             },
    #         ],
    #     )
    #     app = JupyterDash("Network Plot")
    #     app.layout = html.Div([
    #         cy,
    #         html.Div([
    #             html.Div(dcc.Input(id="input-rxn", type="text", placeholder="Enter a reaction"),
    #                 style={"height": 15, "padding": "5px"}),
    #             html.Div(dcc.Input(id="input-const-rxn", type="number", placeholder="Enter a reaction constant"), style={"padding": "5px"}),
    #             html.Button("Add species", id="button-rxn", style={"padding": "5px"}),
    #             html.Div(id='label-rxn', children="", style={"padding": "5px"})
    #         ], style={"display": "flex", "flex-direction": "row"}),
    #         html.Pre(id='cytoscape-tapNodeData-json'),
    #     ])
    #
    #     # Setup a callback in case the edit inputs are used
    #     @app.callback(
    #         dash.dependencies.Output('label-rxn', 'children'),
    #         [dash.dependencies.Input("button-rxn", "n_clicks")],
    #         [dash.dependencies.State("input-rxn", "value")],
    #         [dash.dependencies.State("input-const-rxn", "value")],
    #     )
    #     def add_reaction(n_clicks, rxn, k):
    #         if rxn is None or k is None:
    #             return
    #
    #         # Store the updated positions before adding the new elements.
    #         # TODO(rithvik): This is very brittle
    #         # elems = cy.elements
    #         # for e in elems:
    #             # # Only process vertices, not edges (as those are not changeable through drag-drop).
    #             # if "position" not in e:
    #                 # continue
    #             # d = e["data"]
    #             # p = e["position"]
    #             # print("updating", d["id"], p["x"], p["y"])
    #             # self.network_graph_pos[d["id"]] = np.array([p["x"], p["y"]])
    #         # print(self.network_graph_pos)
    #
    #         success, msg = self._add_raw_reaction(rxn, float(k))
    #         return "Success!" if success else msg
    #
    #     app.run_server(mode="inline")
    #
    # def _add_raw_reaction(self, raw: str, k: float) -> Tuple[bool, str]:
    #     """Add a new reaction given the raw text form and reaction constant.
    #
    #     The species in the reaction must already be defined in the system. True is returned iff the
    #     addition is successful. If the addition fails, a string error message is also returned.
    #
    #     The given raw text is split into reactants and products sections which are individually
    #     parsed to determine matching symbols which can be combined together to create a new
    #     reaction.
    #     """
    #     # Use Species Manager instead of the current class's symbol list as any defined symbol can
    #     # be used in a reaction.
    #     syms = {str(s): s for s in self.species_manager.all_symbols}
    #     sides = [s.strip() for s in raw.split("=")]
    #     if len(sides) != 2:
    #         return False, "Invalid equation format: there must be exactly one \"=\"."
    #
    #     raw_reactants = [s.strip() for s in sides[0].split("+")]
    #     raw_products = [s.strip() for s in sides[1].split("+")]
    #
    #     # Parse a single term such as "4x" into a coefficient and symbol
    #     def parse_term(t: str) -> Tuple[str, int, bool]:
    #         if len(t) < 1:
    #             return "", 0, False
    #         i = 0
    #         c = t[i]
    #         raw_coeff = ""
    #         while c.isdigit():
    #             raw_coeff += c
    #             i += 1
    #             c = t[i]
    #         if raw_coeff == "":
    #             raw_coeff = "1"
    #         return t[i:], int(raw_coeff), True
    #
    #     # Create lists of tuples of sympy symbols and coefficients
    #     reactants = []
    #     for r in raw_reactants:
    #         v, coeff, ok = parse_term(r)
    #         if not ok:
    #             return False, f"Invalid reactant species {r}."
    #         if v not in syms:
    #             return False, f"Reactant species does not exist {v}."
    #         reactants.append((syms[v], coeff))
    #
    #     products = []
    #     for p in raw_products:
    #         v, coeff, ok = parse_term(p)
    #         if not ok:
    #             return False, f"Invalid product species {p}."
    #         if v not in syms:
    #             return False, f"Product species does not exist {v}."
    #         products.append((syms[v], coeff))
    #
    #     if len(reactants) == 0:
    #         return False, "At least one reactant is required"
    #
    #     # Create the sympy expressions
    #     inp = reactants[0][1] * reactants[0][0]
    #     for i in range(1, len(reactants)):
    #         inp += reactants[i][1]*reactants[i][0]
    #
    #     out = products[0][1] * products[0][0]
    #     for i in range(1, len(products)):
    #         inp += products[i][1]*products[i][0]
    #
    #     rxn = Rxn(inp, out, k)
    #
    #     self.reactions.append(rxn)
    #     self.terms.extend(rxn.to_terms())
    #     self._update_symbols()
    #
    #     self._add_to_network_graph(rxn.reactants.free_symbols, rxn.products.free_symbols, rxn.rate_constant)
    #
    #     return True, ""

    def display_ode_expressions(self):
        """TODO"""
        try:
            from IPython.core import display
        except ModuleNotFoundError:
            warnings.warn('BulkRxnSystem.display_ode_expressions requires '
                          'IPython.')
            return

        time = sym.symbols('t')
        diffs = self.get_ode_expressions()
        symbols = self.get_symbols_ordered()

        eq_tuples = list(zip(symbols, diffs))
        eqs = []

        for symbol, deriv in eq_tuples:
            eqs.append(sym.Eq(sym.Derivative(symbol, time), deriv,
                              evaluate=False))
        for eq in eqs:
            display.display(eq)

    @util.depreciate
    def text(self) -> str:
        """Return a text representation of the reaction system, describing the chemical equations in
        natural language."""
        text: str = ""
        for rxn in self.reactions:
            text += rxn.text() + " "
        return text[:-1]

    def __str__(self):
        s = self.__class__.__name__ + ' with components:\n'
        for component in self.elements:
            comp_lines = str(component).splitlines()
            s += ''.join([f'\t{line}\n' for line in comp_lines])
        return s[:-1]

    def __repr__(self):
        return f'{self.__class__.__name__}(components={repr(self.elements)})'

    @util.depreciate
    def id(self):
        """Return a unique identifier for the reactions in this system.

        This identifier does not ignore concentrations or reaction constants. The order of reactions
        (and whether they are formed as reversible or pairs of regular reactions) also do not
        matter.
        """
        f = []
        for c in self.elements:
            if isinstance(c, BulkRxn):
                f.extend(c.id())
            elif isinstance(c, conditions.Schedule):
                f.append(repr(c)) # TODO(rithvik): This is a hack
        return "_".join(sorted(f))

    @util.depreciate
    def fingerprint(self):
        """Return a unique fingerprint for the reactions in this system.

        This identifier ignores concentrations and reaction constants. The order of reactions
        (and whether they are formed as reversible or pairs of regular reactions) also do not
        matter.
        """
        f = []
        for c in self.elements:
            if isinstance(c, BulkRxn):
                f.extend(c.fingerprint())
        return "_".join(sorted(f))

    @property
    @util.depreciate
    def ode(self):
        """
        Return the list of ODE expressions.
        """
        return self.get_ode_expressions()

    @util.depreciate
    def show_ode(self):
        """
        Print the ODEs line by line.
        """
        for e in self.ode:
            print(e)

    @util.depreciate
    def tp_enum(self, dft_outputs: Optional[Dict[str, str]] = None):
        """
        Build a new TP relation class based on DFT calculations.

        :dft_ios: a dictionary from reaction name to the corresponding DFT output file path.
        """
        self.tprate_relation = {}
        for name, dft_io in dft_outputs.items():
            rxn = self.elements_by_name[name]
            if isinstance(rxn, BulkRevRxn):
                self.tprate_relation[rxn.name] = TPRateRelation(tp_file=dft_io,
                                                                init_temp=self.temperature.value,
                                                                init_pressure=self.pressure.value,
                                                                # TODO: bug-proof this
                                                                adsorption=[0] if rxn.is_adsorption(self.species_manager) else [1],
                                                                desorption=[0] if rxn.is_desorption(self.species_manager) else [1],
                                                                # TODO: include Surface Rxns
                                                                init_constants=[rxn.rate_constant, rxn.rate_constant])
            else:
                self.tprate_relation[rxn.name] = TPRateRelation(tp_file=dft_io,
                                                                init_temp=self.temperature.value,
                                                                init_pressure=self.pressure.value,
                                                                adsorption=[0] if rxn.is_adsorption(self.species_manager)  else [],
                                                                desorption=[0] if rxn.is_desorption(self.species_manager)  else [],
                                                                # TODO: include Surface Rxns
                                                                init_constants=[rxn.rate_constant])
        self.tpconc_relation = {}
        for s in self.schedules:
            if self.species_manager.is_gas(s.symbol):
                self.tpconc_relation[s.symbol] = IdealGasLawRelation(p_in_torr=self.pressure.value,
                                                                     n=s.initial_concentration,
                                                                     t=self.temperature.value)

    @util.depreciate
    def tp_next_rsys(self,
                     t=None,
                     total_p=None,
                     partial_pressures=None,
                     include_rules=False,
                     inplace=False):
        """
        Return a new Rsys following same TP relation.


        :param partial_pressures: a dictionary from a gas species to its partial pressure.
                                  if total_p is set, partial_pressure is a percentage of its value;
                                  otherwise, use the original temperature value;
                                  if a gas species does not appear in this dictionary, we assume it has same pressure
                                  as the partial pressure.
        :param include_rules: if set to True, apply partial pressures to adsorptions and desorptions involving the only
                              relevant gas.
        :return: if inplace, return None; otherwise, return a RxnSystem.
        """
        if total_p is None:
            total_p = self.pressure.value
        if t is None:
            t = self.temperature.value

        # Ensure that partial pressure is a dictionary from species name to a floating number.
        new_partial_pressures = defaultdict(lambda: 1)
        if partial_pressures:
            for k, v in partial_pressures.items():
                if isinstance(k, sym.Symbol):
                    new_partial_pressures[k.name] = v
                else:
                    new_partial_pressures[k] = v
        partial_pressures = new_partial_pressures

        # TODO: check this step
        next_constants = {}
        for name, r in self.tprate_relation.items():
            # TODO: try to trouble shoot the following lines.
            sorption_species = self.elements_by_name[name].sorption_species(sm=self.species_manager)
            if include_rules and sorption_species.name in partial_pressures:
                pressure = total_p * partial_pressures[sorption_species.name]
                if pressure == 0:
                    pressure = 1e-15
                print(f"{str(self.elements_by_name[name])} at partial pressure {partial_pressures[sorption_species.name]}")
            else:
                pressure = total_p


            next_constants[name] = r.constants(temp=t, pressure=pressure)
        next_rsys = self if inplace else copy.deepcopy(self)
        for name, next_constants in next_constants.items():
            rxn = next_rsys.elements_by_name[name]
            if rxn.is_reversible:
                rxn.set_rates(rate=next_constants[0], rate_reverse=next_constants[1])
            else:
                rxn.set_rate(rate=next_constants[0])
        next_rsys.temperature.value = t
        next_rsys.pressure.value = total_p

        # Reset gas molecular count, assuming that volume is constant.
        for s in next_rsys.schedules:
            if next_rsys.species_manager.is_gas(s.symbol):
                partial_pressure = total_p * partial_pressures[s.symbol.name]
                print(f"{s.symbol} has a partial pressure of {partial_pressures[s.symbol.name]}")


                new_n = next_rsys.tpconc_relation[s.symbol].calculate_n(p=partial_pressure, t=t)
                s.update_conc(func=lambda x: new_n, inplace=True)
        if inplace:
            return
        else:
            return next_rsys

