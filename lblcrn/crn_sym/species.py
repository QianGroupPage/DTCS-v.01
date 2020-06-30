"""
CRN - species.py:

Structures for the manipulation of chemical speices.

Credits:
Dr. Jin Qian, Domas Buracas, Ye Wang, Andrew Bogdan, Rithvik Panchapakesan
"""

# *** Libraries ***
import sympy as sym
# TODO: one way to enforce order in the expressions is to overwrite add's sorting function to
#  preversve order. However, this is highly dependent on sympy.
# sym.core.Add._addsort = lambda x: None
from lblcrn.common import color_to_RGB
from lblcrn.crn_sym.surface import Site

from typing import List, Tuple, Union

# *** Classes ***
class Orbital:
    """
    An orbital in a species, this is essentially a named tuple, it's a class for readability purposes.
    """

    def __init__(self, name: str, binding_energy: float, splitting: float=1):
        self.name = name
        self.binding_energy = binding_energy
        self.splitting = splitting

    def __str__(self):
        if self.splitting == 1:
            return self.name + '@ ' + str(self.binding_energy)
        else:
            return self.name + '@ ' + repr(self.binding_energy) + ', splitting ' + repr(self.splitting)

    def __repr__(self):
        if self.splitting == 1:
            return "Orbital(name=" + self.name + ', binding_energy=' + str(self.binding_energy) + ')'
        else:
            return "Orbital(name=" + self.name + ', binding_energy=' + repr(self.binding_energy) + ', splitting=' + repr(self.splitting) + ')'


class Species:
    """
    A chemical species with a name and orbitals, which are triples of (orbital name, binding energy, proportion)
    """

    def __init__(self, name: str, orbitals: List[Orbital], color: Union[Tuple[int], List[int], str] = None,
                 parent=None, site: Site = None):
        self.name = name
        self.orbitals = orbitals
        if color:
            self.color = color_to_RGB(color)
        else:
            self.color = color
        self.parent = parent   # The parent of the species
        self.sub_species = {}   # The dictionary of sub species from name to the object

        self.site = site
        if self.site and self.site != Site.default:
            self.create_sub_species(suffix=site.name, color=self.color)

    def create_sub_species(self, suffix: str = "", color: Union[Tuple[int], List[int], str] = "", entire_name: str =""):
        if not entire_name:
            if not suffix:
                suffix = f"sub_{len(self.sub_species)}"
            entire_name = f"{self.name}_{suffix}"
        elif suffix:
            raise Exception(f"Both suffix={suffix} and entire_name {entire_name} provided to " +
                            f"create_sub_species_function")
        if not color:
            color = self.color
        sub_species = Species(entire_name, self.orbitals, color=color, parent=self)

        if entire_name in self.sub_species:
            raise Exception(f"species {self.name} already has sub species {repr(self.sub_species[entire_name])} " +
                            f"with name {entire_name}.")
        else:
            self.sub_species[entire_name] = sub_species
        return

    def __str__(self):
        if self.color:
            return self.name + ", orbitals: " + str(self.orbitals) + ", color: " + str(self.color)
        return self.name + ", orbitals: " + str(self.orbitals)

    def __repr__(self):
        return 'Species(name=' + self.name + ', orbitals=' + repr(self.orbitals) + ', size=' + repr(self.color) + ')'


class SpeciesManager:
    """
    A smart wrapper of a dictionary {sym.Symbol: Species} for the purpose of keeping track of
    which symbols correspond to which speices.

    You can create symbols/species pairs with SpeciesManager.sp and access them with SpeciesManager[],
    which forward to the more verbosely-named make_species and species_from_symbol
    """

    def __init__(self):
        self._species = {} # As of current, initializes empty

    def make_species(self, name: str, orbitals: Union[Orbital, List[Orbital]], site: Site = None) -> sym.Symbol:
        """
        Makes a sym.Symbol and a corresponding Species and keeps track of their correspondence.
        Returns the symbol.

        Orbitals can be either a list of orbitals or just one orbital, just to be kind.
        """
        symbol = sym.Symbol(name)

        if not isinstance(orbitals, list):
            orbitals = [orbitals]

        self._species[symbol] = Species(name, orbitals, site=site)
        return symbol

    def species_from_symbol(self, key: sym.Symbol) -> Species:
        return self._species[key]

    def __str__(self):
        return str(self._species)  # TODO

    def __repr__(self):
        pass  # TODO

    sp = make_species
    __getitem__ = species_from_symbol
