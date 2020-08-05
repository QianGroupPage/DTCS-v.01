"""
CRN - surface.py:

Structures for the creation of surfaces

Credits:
Dr. Jin Qian, Domas Buracas, Andrew Bogdan, Rithvik Panchapakesan, Ye Wang
"""

# *** Libraries ***
from typing import List, Tuple, Union
import sympy as sym
from lblcrn.common import color_to_RGB
from lblcrn.surface_crn.surface_crns.models.coord_grid import CoordGrid

# *** Classes ***
class Surface:
    """
    A surface structure, by default, it is square with only top sites.
    """

    def __init__(self, name: str, size: Tuple[int] = (0, 0), structure: str = "rectangle",
                 color: Union[Tuple[int], List[int], str] = None, poscar_file: str = "", supercell_dimensions=1):
        # This surface is based on a POSCAR file.
        if poscar_file:
            self.use_coord_grid = True
            self.poscar_file = poscar_file
            self.supercell_dimensions = supercell_dimensions
            self.coord_grid = CoordGrid.from_poscar(poscar_file, supercell_dimensions=supercell_dimensions)

            # a dummy variable for coord_grid's 2-d size
            size = self.coord_grid.size_2d
        else:
            self.use_coord_grid = False

        self.name = name
        self.size = size
        self.structure = structure

        if color:
            self.color = color_to_RGB(color)
        else:
            self.color = color

        if self.structure not in Surface.allowed_structures():
            raise Exception(f"Structure {structure} is not allowed.\n\n" +
                            f"The only allowed sites are " +
                            " ".join(Surface.allowed_structures()))

        self.populate_sites()  # Populate all the sites for this surface

    @staticmethod
    def allowed_structures():
        return "rectangle", "hexagon"

    @property
    def allowed_sites(self):
        if self.use_coord_grid:
            # TODO: build "fold names"
            return self.coord_grid.fold_names
        elif self.structure == "rectangle":
            return ["top"]
        elif self.structure == "hexagon":
            return ["top", "threefold"]

    @property
    def sites(self):
        res = []
        for s in self.allowed_sites:
            if s != "top":
                res.append(getattr(self, s))
        return res

    @property
    def top(self):
        """
        A surface itself can represent its top site.
        :return: the object itself
        """
        return self

    def populate_sites(self):
        for site_name in self.allowed_sites:
            if site_name != "top":
                setattr(self, site_name, Site(site_name, self))

            # TODO: verify the correct modifications for 3F sites.
            # if site_name == "threefold":
            #     setattr(self, site_name, Site("3F", self))
            # elif site_name != "top":
            #     setattr(self, site_name, Site(site_name, self))

    # @property
    def symbol(self):
        return sym.Symbol(self.name)

    @property
    def symbols(self):
        symbols = set([self.symbol()])
        for s in self.sites:
            symbols.add(sym.Symbol(s.name))
        return symbols

    def __str__(self):
        row = " ".join([self.name] * self.size[1])
        if self.color:
            return "\n".join([row] * self.size[0]) + f"\n color={self.color}\n" 
        return "\n".join([row] * self.size[0]) 

    def __repr__(self):
        return "Surface(name=" + self.name + ', size=' + repr(self.size) + " , color=" + repr(self.color) + ')'


class Site:
    """
    A site on a surface
    """
    default = "top"

    def __init__(self, name: str, surface: Surface, color: Union[Tuple[int], List[int], str] = None):
        self.name = name
        self.surface = surface

        if color:
            self.color = color_to_RGB(color)
        else:
            self.color = color

    @property
    def symbol(self):
        return sym.Symbol(self.name)

    def set_color(self, color):
        self.color = color_to_RGB(color)

