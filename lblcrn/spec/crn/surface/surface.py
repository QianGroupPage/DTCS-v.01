"""
CRN - surface.py:

Structures for the creation of surfaces

Credits:
Dr. Jin Qian, Domas Buracas, Andrew Bogdan, Rithvik Panchapakesan, Ye Wang
"""

# *** Libraries ***
from __future__ import annotations

from typing import List, Tuple, Union

import random

import sympy as sym

from lblcrn.common import color_to_RGB, util
from lblcrn.common.num_to_word import num2word
from lblcrn.sim.surface_crn.connectivity.triangulation import grid_size, show_triangulation
from lblcrn.sim.surface_crn.connectivity.voronoi import fold_numbers, produce_voronoi
from lblcrn.spec.spec_abc import Spec
from lblcrn.common.colors import color_map

from lblcrn.sim.surface_crn.surface_crns.models.grids import *
from lblcrn.sim.surface_crn.surface_crns.views.grid_display import *
# from lblcrn.sim.surface_crn.surface_crns.models.grids import SquareGrid
# from lblcrn.sim.surface_crn.api_adapter.api_adapt import HexGridPlusIntersections


# *** Classes ***
class Surface(Spec):
    """
    A surface structure, by default, it is square with only top sites.
    """

    valid_structures = ['rectangle', 'hexagon']

    def __init__(self,
                 name: str,
                 # size: Tuple[int] = (10, 10), TODO(Andrew) you _should_ be able to supply this, just preffered later
                 structure: str = 'rectangle',
                 color: Color = None,
                 poscar_file: str = '',
                 # supercell_dimensions=1,
                 # surface_depth=1,
                 description=''):
        self.name = name
        # self.size = size
        self.structure = structure

        # A dictionary from site names in allowed_sites to nested dictionaries from species names to the
        # initial concentration of the species.
        self.initial_species = None

        if color:
            color_map[self.name] = color
            self.color = color

        if self.structure not in self.valid_structures:
            raise Exception(f'Structure {structure} is not allowed.\n\n' +
                            f'The only allowed sites are ' +
                            ' '.join(self.valid_structures))

        self._make_sites()

        # self.populate_sites()  # Populate all the sites for this surface
        super().__init__(name)

    @classmethod
    def from_poscar(cls, poscar):
        raise NotImplementedError()
        # This surface is based on a POSCAR file.
        self.structure = 'poscar'
        if poscar:
            self.use_coord_grid = True
            self.poscar = poscar
            self.supercell_dimensions = supercell_dimensions
            self.surface_depth = surface_depth
            self.default_names = {}

            _, points = show_triangulation(poscar_file)
            self.vor = produce_voronoi(points)

            # a dummy variable for coord_grid's 2-d size
            size = grid_size(points)

    # def set_initial_concentrations(self, species):
    #     self.initial_species = species
    #
    # def set_default_name(self, site_name, species_name):
    #     self.default_names[site_name] = species_name

    def _make_sites(self):
        if self.structure == 'poscar':
            # TODO: do we always have two fold sites?
            sites = ['top', 'twofold'] + [f'{num2word(num)}fold' for num in fold_numbers(self.vor)]
            # return self.coord_grid.fold_names
        elif self.structure == 'rectangle':
            sites = ['top']
        elif self.structure == 'hexagon':
            sites = ['top', 'threefold']

        self.sites = []
        for site in sites:
            site_name = f'{self.name}_{site}'
            self.sites.append(site_name)
            setattr(self, site, sym.Symbol(site_name))
            color_map[site_name] = color_map[self.name]  # TODO(Andrew): Give more options for this

    def make_state(self, *coverages, size=(10, 10)):
        if self.structure == 'rectangle':
            return self._make_state_rectangle(coverages, size)
        elif self.structure == 'hexagon':
            return self._make_state_hexagon(coverages, size)
        assert False, 'Unreachable'

    def _make_state_rectangle(self, coverages, size):
        assert hasattr(self, 'top'), 'Rectangle should have a top site'

        rows = size[0]
        cols = size[1]

        state = [[self.top.name for _ in range(cols)] for _ in range(rows)]
        for row in range(rows):
            for col in range(cols):
                for coverage in coverages:
                    # TODO(Andrew): this doesn't actually give correct coverage, due to over-writing
                    if coverage.coverage < random.random():
                        state[row][col] = coverage.species

        # grid = SquareGrid(rows, cols)
        # grid.set_global_state(state)

        return np.asarray(state)

    def _make_state_hexagon(self, coverages, size):
        assert hasattr(self, 'top'), 'Hex should have a top site'
        assert hasattr(self, 'threefold'), 'Hex should have a threefold site'

        rows = size[0]
        cols = size[1]

        state = [[self.top.name for _ in range(cols)] for _ in range(rows)]
        for row in range(rows):
            for col in range(cols):
                for coverage in coverages:
                    if coverage.site != self.top.name:
                        continue
                    # TODO(Andrew): this doesn't actually give correct coverage, due to over-writing
                    if coverage.coverage < random.random():
                        state[row][col] = coverage.species

        state_3f = [[self.threefold.name for _ in range(rows - 1)]
                    for _ in range((rows - 1) * 2)]
        for row in range(rows):
            for col in range(cols):
                for coverage in coverages:
                    if coverage.site != self.threefold.name:
                        continue
                    # TODO(Andrew): this doesn't actually give correct coverage, due to over-writing
                    if coverage.coverage < random.random():
                        state[row][col] = coverage.species

        return np.asarray(state), np.asarray(state_3f)

        # species = {"threefold": [], "top": []}
        # for c in rsys.schedules:
        #     # TODO
        #     if isinstance(c, Schedule) and not isinstance(c, Conc):
        #         raise Exception(f"{c} is a schedule, not an initial concentration." +
        #                         f" This is not currently supported in the Surface CRN.")
        #
        #     if not rsys.species_manager.species_from_symbol(c.symbol).site or \
        #             rsys.species_manager.species_from_symbol(c.symbol).site.name == "top":
        #         species["top"].extend([str(c.symbol) for _ in range(int(str(c.concentration)))])
        #     elif rsys.species_manager.species_from_symbol(c.symbol).site.name == "threefold":
        #         species["threefold"].extend([str(c.symbol) for _ in range(int(str(c.concentration)))])
        #
        # # Initiate the hexagonal grid
        # rows = rsys.surface.size[0]
        # cols = rsys.surface.size[1]
        # surface = HexGridPlusIntersections(rows, cols)
        # for n in surface:
        #     if n.is_intersection:
        #         n.state = "threefold"
        #     else:
        #         n.state = rsys.surface.name
        #
        # # Popultate the structure with initial atoms.
        # initial_nodes = {"threefold": [], "top": []}
        # for n in surface:
        #     if n.is_intersection:
        #         initial_nodes["threefold"].append(n)
        #     else:
        #         initial_nodes["top"].append(n)
        #
        # random.seed(random_seed)
        # for k, v in species.items():
        #     nodes = initial_nodes[k]
        #     indices = random.sample(range(len(nodes)), len(v))
        #     chosen = [nodes[i] for i in indices]
        #     for i, n in enumerate(chosen):
        #         n.state = species[k][i]
        #
        # return surface

    # @property
    # def sites(self):
    #     res = []
    #     for s in self.allowed_sites:
    #         if s != 'top':
    #             res.append(getattr(self, s))
    #     return res

    # @property
    # def top(self):
    #     """
    #     A surface itself can represent its top site.
    #     :return: the object itself
    #     """
    #     return self

    # def populate_sites(self):
    #     for site_name in self.allowed_sites:
    #         if site_name != 'top':
    #             setattr(self, site_name, Site(site_name, self))
    #
    #         # TODO: verify the correct modifications for 3F sites.
    #         # if site_name == 'threefold':
    #         #     setattr(self, site_name, Site('3F', self))
    #         # elif site_name != 'top':
    #         #     setattr(self, site_name, Site(site_name, self))

    # # @property
    # def symbol(self):
    #     return sym.Symbol(self.name)

    # @property
    # def symbols(self):
    #     symbols = set([self.symbol()])
    #     for s in self.sites:
    #         symbols.add(sym.Symbol(s.name))
    #     return symbols

    # def __str__(self):
    #     row = ' '.join([self.name] * self.size[1])
    #     if self.color:
    #         return '\n'.join([row] * self.size[0]) + f'\n color={self.color}\n'
    #     return '\n'.join([row] * self.size[0])
    #
    # def __repr__(self):
    #     return 'Surface(name=' + self.name + ', size=' + repr(self.size) + ' , color=' + repr(self.color) + ')'


class Site:
    """
    A site on a surface
    """
    default = 'top'

    @util.depreciate
    def __init__(self, name: str, surface: Surface, color: Union[Tuple[int], List[int], str] = None):
        self.name = name
        self.surface = surface

        color_map[self.name] = color
        if color:
            self.color = color_to_RGB(color)
        else:
            self.color = color

    @property
    def symbol(self):
        return sym.Symbol(self.name)

    def set_color(self, color):
        self.color = color_to_RGB(color)
