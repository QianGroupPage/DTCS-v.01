import copy
import random
from typing import List, Tuple, Dict

import monty.json
import networkx as nx
import sympy as sym
import plotly.graph_objects as go
import numpy as np

from jupyter_dash import JupyterDash
import dash_cytoscape as cyto
import dash_html_components as html
import dash_core_components as dcc
import dash
import IPython.display as display

import lblcrn
from lblcrn.common import color_to_RGB, generate_new_color
from lblcrn.crn_sym import conditions, species, surface
from lblcrn.crn_sym.reaction import RevRxn, Rxn
from lblcrn.crn_sym.surface_reaction import SurfaceRxn


class RxnSystem(monty.json.MSONable):
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

    def __init__(self, *components):
        """Create a new reaction system. Requires a SpeciesManager.

        Accepts Rxns, Revrxns, Concs, Schedules, Terms, ConcEqs,
        ConcDiffEqs, and (one) SpeciesManager in any order.

        If you have a function returning a collection of the above, you do
        not have to worry about unpacking the collection: it will unpack and
        flatten lists and tuples for you.
        """
        assert len(components) > 0, 'Must pass at least one reaction and a species manager to a solution system.'
        if len(components) == 1:
            assert not isinstance(components[0], species.SpeciesManager), 'Must pass at least one reaction to a solution system.'
            raise AssertionError('Must pass a species manager to a solution system.')

        # Flatten the components
        flat = False
        while not flat:
            flatter_components = []
            flat = True
            for component in components:
                if isinstance(component, list) or isinstance(component, tuple):
                    flat = False
                    flatter_components.extend(component)
                else:
                    flatter_components.append(component)
            components = flatter_components
        self.components = flatter_components  # pyint: disable=

        # Split into terms, schedules, and conc (diff.) eq.s
        # These are not sorted by the symbol index.
        self.terms = [] #: Terms are terms.
        self.surface_rxns = []
        self.schedules = []
        self.conc_eqs = []
        self.conc_diffeqs = []
        self.species_manager = None
        self.surface = None
        self.reactions = []

        for component in self.components:
            if isinstance(component, conditions.Schedule):
                self.schedules.append(component)
            elif isinstance(component, SurfaceRxn):
                # self.terms.extend(component.to_terms())
                # pass
                self.surface_rxns.append(component)
            elif isinstance(component, Rxn):
                self.reactions.append(component)
                self.terms.extend(component.to_terms())
            elif isinstance(component, conditions.Term):
                self.terms.append(component)
            elif isinstance(component, conditions.ConcEq):
                self.conc_eqs.append(component)
            elif isinstance(component, conditions.ConcDiffEq):
                self.conc_diffeqs.append(component)
            elif isinstance(component, surface.Surface):
                self.surface = component
            elif isinstance(component, species.SpeciesManager):
                self.species_manager = component
            else:
                raise AssertionError(f'Unknown input {component} of type ' + str(type(component)) + ' to reaction system.')

        assert self.species_manager is not None, 'Must pass a species manager to a solution system.'

         # Share the surface name to the species manager
        if self.surface is not None:
            self.species_manager.default_surface_name = self.surface.name

        self._update_symbols()
        self._generate_network_graph()

    def _update_symbols(self):
        """Iterate over reactions etc and update the symbols set and scheduler.

        This method should be called whenever the reaction system is updated.
        """
        # Pick an order for the symbol
        self._symbols = set()
        for rxn in self.surface_rxns:
            for s in rxn.get_symbols():
                if s not in self.surface.symbols:
                    self._symbols.add(s)
        for term in self.terms:
            self._symbols.update(term.get_symbols())
        for schedule in self.schedules:
            self._symbols.add(schedule.symbol)
        for equation in self.conc_eqs:
            self._symbols.update(equation.get_symbols())
        for equation in self.conc_diffeqs:
            self._symbols.update(equation.get_symbols())
        self._symbols = sorted(list(self._symbols), key=lambda s: str(s))

        # Make an indexing dictionary
        self.symbol_index = {}
        for index, symbol in enumerate(self._symbols):
            self.symbol_index[symbol] = index

        # Make symbol:concentration/schedule list
        # Make default (Conc 0 @ t=0 for each species) scheduler list
        self.scheduler = [conditions.Conc(symbol, 0) for symbol in self._symbols]
        # Overwrite scheduler with Concs and Schedules
        for schedule in self.schedules:
            self.scheduler[self.symbol_index[schedule.symbol]] = schedule

        # an index from symbols to their colors
        self.color_index = None

    def get_ode_expressions(self) -> List[sym.Expr]:
        """Return a list of expressions, corresponding to the derivative of the
        concentration of each symbol in the reaction system.

        The collection is ordered, and that order is accessible through
        symbol_index or through get_symbols.
        """

        # Make an emtpy ODE list
        odes = [sym.sympify(0)] * len(self._symbols)

        # Sum all terms for each symbol
        for term in self.terms:
            odes[self.symbol_index[term.symbol]] += term.expression

        # Set the conc diffeqs
        for equation in self.conc_diffeqs:
            odes[self.symbol_index[equation.symbol]] = equation.expression

        return odes

    def get_ode_functions(self):
        """Return the ODE function with signature func(t, y) of the system."""

        symbols = self.get_symbols()
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

    def get_conc_functions(self):
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

    def _get_species(self) -> List[species.Species]:
        species = []

        for symbol in self._symbols:
            species.append(self.species_manager[symbol])

        return species

    def get_symbols(self) -> List[sym.Symbol]:
        """
        Give all the symbols in the reaction system in a fixed order.
        """
        return copy.copy(self._symbols)

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = {
            '@module': self.__class__.__module__,
            '@class': self.__class__.__name__,
            '@version': lblcrn.__version__,  # TODO: Better way to do this?
            'components': [comp.as_dict() for comp in self.components]
        }
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        decode = monty.json.MontyDecoder().process_decoded
        components = [decode(comp) for comp in d['components']]
        return cls(*components)

    # TODO
    def get_colors(self):
        """Assign (if applicable) and return colors for all species.
        """
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

    @property
    def surface_names(self) -> List[str]:
        """A list for names that appear on the surface.
        """
        return [self.surface.name] + [s.name for s in self.surface.sites]

    def _add_to_network_graph(self, r, p, w):
        """Add the given reactions and products to the network graph.
        """
        for reactant in r:
            for product in p:
                self.network_graph.add_edge(str(product), str(reactant), weight=w)

    def _generate_network_graph(self):
        """Create a reaction network graph (data structure).

        Each species becomes a vertex and reactions are represented with directed edges from
        reactants to products (where the reaction constants are the edge weights).
        """
        self.network_graph = nx.DiGraph()

        for rxn in self.reactions:
            r = rxn.reactants.free_symbols
            p = rxn.products.free_symbols
            
            self._add_to_network_graph(r, p, rxn.rate_constant)
            if isinstance(rxn, RevRxn):
                self._add_to_network_graph(p, r, rxn.rate_constant_reverse)
        self.network_graph_pos = nx.spring_layout(self.network_graph)
        for k, pos in self.network_graph_pos.items():
            # Scale positions to make graphing easier.
            self.network_graph_pos[k] = np.array([500*pos[0], 500*pos[1]])
    
    def plot(self):
        """Plot the reaction network as an interactive graph.

        The reaction network graph is plotted in a user-draggable view. Users can also add new
        reactions to the system using the provided text input (simply re-run the cell to view the
        updated graph with the new reaction).
        """
        G = self.network_graph
        pos = self.network_graph_pos
        largest_weight = -1

        # Convert the networkx graph to a format that Cytoscape can use
        for node in G.nodes():
            for e in G.out_edges(node):
                weight = G.get_edge_data(node, e[1])["weight"]
                largest_weight = max(largest_weight, weight)

        elements = []
        for node in G.nodes():
            x, y = pos[node]
            elements.append({"data": {"id": str(node), "label": str(node)}, "position": {"x": x, "y": y}})
            for e in G.out_edges(node):
                weight = G.get_edge_data(node, e[1])["weight"]
                elements.append({"data": {
                    "source": str(node), "target": str(e[1]), "weight": weight, "normalized_weight": weight / largest_weight * 3.75 + 2
                }})

        # Display the plot as well as the related edit inputs
        cy = cyto.Cytoscape(
            id="network-plot",
            layout={"name": "preset"},
            style={"width": "100%", "height": "500px"},
            elements=elements,
            stylesheet=[
                {
                    "selector": "edge",
                    "style": {
                        # The default curve style does not work with certain arrows
                        "curve-style": "bezier",
                        "source-arrow-color": "black",
                        "source-arrow-shape": "triangle",
                        "label": "data(weight)",
                        "width": "data(normalized_weight)",
                    }
                },
                {
                    "selector": "node",
                    "style": {
                        "label": "data(label)",
                    }
                },
            ],
        )
        app = JupyterDash("Network Plot")
        app.layout = html.Div([
            cy,
            html.Div([
                html.Div(dcc.Input(id="input-rxn", type="text", placeholder="Enter a reaction"),
                    style={"height": 15, "padding": "5px"}),
                html.Div(dcc.Input(id="input-const-rxn", type="number", placeholder="Enter a reaction constant"), style={"padding": "5px"}),
                html.Button("Add species", id="button-rxn", style={"padding": "5px"}),
                html.Div(id='label-rxn', children="", style={"padding": "5px"})
            ], style={"display": "flex", "flex-direction": "row"}),
            html.Pre(id='cytoscape-tapNodeData-json'),
        ])

        # Setup a callback in case the edit inputs are used
        @app.callback(
            dash.dependencies.Output('label-rxn', 'children'),
            [dash.dependencies.Input("button-rxn", "n_clicks")],
            [dash.dependencies.State("input-rxn", "value")],
            [dash.dependencies.State("input-const-rxn", "value")],
        )
        def add_reaction(n_clicks, rxn, k):
            if rxn is None or k is None:
                return
            
            # Store the updated positions before adding the new elements.
            # TODO(rithvik): This is very brittle
            # elems = cy.elements
            # for e in elems:
                # # Only process vertices, not edges (as those are not changeable through drag-drop).
                # if "position" not in e:
                    # continue
                # d = e["data"]
                # p = e["position"]
                # print("updating", d["id"], p["x"], p["y"])
                # self.network_graph_pos[d["id"]] = np.array([p["x"], p["y"]])
            # print(self.network_graph_pos)

            success, msg = self._add_raw_reaction(rxn, float(k))
            return "Success!" if success else msg
        
        app.run_server(mode="inline")

    def _add_raw_reaction(self, raw: str, k: float) -> Tuple[bool, str]:
        """Add a new reaction given the raw text form and reaction constant.
        
        The species in the reaction must already be defined in the system. True is returned iff the
        addition is successful. If the addition fails, a string error message is also returned.

        The given raw text is split into reactants and products sections which are individually
        parsed to determine matching symbols which can be combined together to create a new
        reaction.
        """
        # Use Species Manager instead of the current class's symbol list as any defined symbol can
        # be used in a reaction.
        syms = {str(s): s for s in self.species_manager.all_symbols}
        sides = [s.strip() for s in raw.split("=")]
        if len(sides) != 2:
            return False, "Invalid equation format: there must be exactly one \"=\"."

        raw_reactants = [s.strip() for s in sides[0].split("+")]
        raw_products = [s.strip() for s in sides[1].split("+")]

        # Parse a single term such as "4x" into a coefficient and symbol
        def parse_term(t: str) -> Tuple[str, int, bool]:
            if len(t) < 1:
                return "", 0, False
            i = 0
            c = t[i]
            raw_coeff = ""
            while c.isdigit():
                raw_coeff += c
                i += 1
                c = t[i]
            if raw_coeff == "":
                raw_coeff = "1"
            return t[i:], int(raw_coeff), True

        # Create lists of tuples of sympy symbols and coefficients
        reactants = []
        for r in raw_reactants:
            v, coeff, ok = parse_term(r)
            if not ok:
                return False, f"Invalid reactant species {r}."
            if v not in syms:
                return False, f"Reactant species does not exist {v}."
            reactants.append((syms[v], coeff))

        products = []
        for p in raw_products:
            v, coeff, ok = parse_term(p)
            if not ok:
                return False, f"Invalid product species {p}."
            if v not in syms:
                return False, f"Product species does not exist {v}."
            products.append((syms[v], coeff))

        if len(reactants) == 0:
            return False, "At least one reactant is required"

        # Create the sympy expressions
        inp = reactants[0][1] * reactants[0][0]
        for i in range(1, len(reactants)):
            inp += reactants[i][1]*reactants[i][0]

        out = products[0][1] * products[0][0]
        for i in range(1, len(products)):
            inp += products[i][1]*products[i][0]

        rxn = Rxn(inp, out, k)

        self.reactions.append(rxn)
        self.terms.extend(rxn.to_terms())
        self._update_symbols()

        self._add_to_network_graph(rxn.reactants.free_symbols, rxn.products.free_symbols, rxn.rate_constant)
        
        return True, ""

    def text(self) -> str:
        """Return a text representation of the reaction system, describing the chemical equations in natural language."""
        text: str = ""
        for rxn in self.reactions:
            text += rxn.text() + " "
        return text[:-1]

    def __str__(self):
        s = self.__class__.__name__ + ' with components:\n'
        for component in self.components:
            comp_lines = str(component).splitlines()
            s += ''.join([f'\t{line}\n' for line in comp_lines])
        return s[:-1]

    def __repr__(self):
        return f'{self.__class__.__name__}(components={repr(self.components)})'
