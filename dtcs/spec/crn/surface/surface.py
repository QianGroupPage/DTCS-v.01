"""
CRN - surface.py:

Structures for the creation of surfaces

Credits:
Dr. Jin Qian, Domas Buracas, Andrew Bogdan, Rithvik Panchapakesan, Ye Wang
"""

# *** Libraries ***
from __future__ import annotations

from typing import List, Tuple, Union, Dict
from dtcs.common.type_alias import SpeciesName, SiteName

import random
import itertools

import sympy as sym

from dtcs.common import util
from dtcs.common.colors.color_gradient import color_to_RGB
from dtcs.common.num_to_word import num2word
from dtcs.sim.surface_crn.connectivity.triangulation import grid_size, show_triangulation
from dtcs.sim.surface_crn.connectivity.voronoi import fold_numbers, produce_voronoi
from dtcs.spec.spec_abc import Spec
from dtcs.common.colors import color_map

from dtcs.sim.surface_crn.surface_crns.models.grids import *
from dtcs.sim.surface_crn.surface_crns.views.grid_display import *
# from dtcs.sim.surface_crn.surface_crns.models.grids import SquareGrid
# from dtcs.sim.surface_crn.api_adapter.api_adapt import HexGridPlusIntersections


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

    def make_state(self,
                   coverage_info: Dict[SiteName, (SpeciesName, float)],
                   size: Tuple[int, int] = (10, 10)):
        if self.structure == 'rectangle':
            return self._make_state_rectangle(coverage_info, size)
        elif self.structure == 'hexagon':
            return self._make_state_hexagon(coverage_info, size)
        assert False, 'Unreachable'

    def _make_state_rectangle(self, coverage_info, size):
        assert hasattr(self, 'top'), 'Rectangle should have a top site'

        rows, cols = size
        sites = list(itertools.product(range(rows), range(cols)))

        sites_map = self._cover_sites(
            sites=sites,
            default=self.top.name,
            coverages=coverage_info[self.top.name]
        )

        state = np.empty(size, np.dtype(object))
        for (row, col), specie in sites_map.items():
            state[row, col] = specie

        return state

    def _make_state_hexagon(self, coverage_info, size):
        assert hasattr(self, 'top'), 'Hex should have a top site'
        assert hasattr(self, 'threefold'), 'Hex should have a threefold site'

        rows_top, cols_top = size
        rows_3f, cols_3f = (rows_top - 1) * 2, cols_top - 1

        sites_top = list(itertools.product(range(rows_top), range(cols_top)))
        sites_3f = list(itertools.product(range(rows_3f), range(cols_3f)))

        sites_map_top = self._cover_sites(
            sites=sites_top,
            default=self.top.name,
            coverages=coverage_info[self.top.name]
        )
        sites_map_3f = self._cover_sites(
            sites=sites_3f,
            default=self.threefold.name,
            coverages=coverage_info[self.threefold.name]
        )

        state_top = np.empty(size, np.dtype(object))
        state_3f = np.empty((rows_3f, cols_3f), np.dtype(object))

        for (row, col), specie in sites_map_top.items():
            state_top[row, col] = specie
        for (row, col), specie in sites_map_3f.items():
            state_3f[row, col] = specie

        return state_top, state_3f

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

    @staticmethod
    def _cover_sites(sites: List,
                     default: SiteName,
                     coverages: List[(SpeciesName, float)]):
        """Cover the given sites according to coverages.

        Args:
            sites: Iterable of any hashable type
            default: str, the un-covered site
            coverages: list of (species, % coverage)

        Returns:
            A dictionary of site: species, where site is pulled from sites,
            and the proportion of each species should roughly match coverages.
        """
        # TODO(Andrew): Raise an error if you exceed 100%?

        sites = list(sites)
        num_sites = len(sites)
        random.shuffle(sites)

        sites_map = {}
        lower = 0
        upper = 0
        for species, coverage in coverages:
            upper = (coverage * num_sites) + lower
            for index in range(math.ceil(lower), math.ceil(upper)):
                sites_map[sites[index]] = species
            lower = upper

        for index in range(math.ceil(upper), num_sites):
            sites_map[sites[index]] = default

        return sites_map

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
