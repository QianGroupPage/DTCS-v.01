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
                 parent=None, site: Site = None, include_sub_species: bool=True):
        self.name = name
        self.orbitals = orbitals
        if color:
            self.color = color_to_RGB(color)
        else:
            self.color = color
        self.parent = parent   # The parent of the species
        self.sub_species = {}   # The dictionary of sub species from name to the object

        self.site = site
        if include_sub_species and self.site and self.site != Site.default:
            self.create_sub_species(suffix=site.name, color=self.color)

    def create_sub_species(self, suffix: str = "", color: Union[Tuple[int], List[int], str] = "", entire_name: str ="",
                           orbitals: Union[Orbital, List[Orbital]] = None, site: Site = None):
        if not entire_name:
            if site:
                suffix = site.name
            elif not suffix:
                suffix = f"sub_{len(self.sub_species)}"
            entire_name = f"{self.name}_{suffix}"
        elif suffix:
            raise Exception(f"Both suffix={suffix} and entire_name {entire_name} provided to " +
                            f"create_sub_species_function")
        if not color:
            color = self.color
        if orbitals is None:
            orbitals = self.orbitals

        sub_species = Species(entire_name, orbitals, color=color, parent=self, site=site, include_sub_species=False)

        if entire_name in self.sub_species:
            raise Exception(f"species {self.name} already has sub species {repr(self.sub_species[entire_name])} " +
                            f"with name {entire_name}.")
        else:
            self.sub_species[entire_name] = sub_species
        return sub_species

    def sub_species_by_name(self, name: str):
        """
        :param name: name of the species
        :return: the species, if found
        """
        if name in self.sub_species:
            return self.sub_species[name]
        else:
            raise Exception(f"Species {self.name} doesn't have a sub-species with name {name}.")

    def sub_species_by_suffix(self, suffix: str):
        """
        :param suffix: the suffix of the species
        :return: the sub_species with the suffix, if found
        """
        return self.sub_species_by_name(self.sub_species_name(suffix))

    def sub_species_name(self, suffix: str):
        """
        :param suffix: suffix of the sub_species
        :return: the sub_species's name
        """
        return f"{self.name}_{suffix}"

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

    def make_species(self, name: Union[str, sym.Symbol], orbitals: Union[Orbital, List[Orbital]] = None,
                     site: Site = None, sub_species_name: str = "") -> sym.Symbol:
        """
        Makes a sym.Symbol and a corresponding Species and keeps track of their correspondence.
        Returns the symbol.

        Orbitals can be either a list of orbitals or just one orbital, just to be kind.
        """
        if isinstance(name, sym.Symbol):
            parent = self.species_from_symbol(name)

            if sub_species_name:
                name = sub_species_name
                if site:
                    name += site.name
                s = parent.create_sub_species(entire_name=name, site=site)
                symbol = sym.Symbol(name)
            else:
                s = parent.create_sub_species(site=site)
                symbol = sym.Symbol(parent.sub_species_name(site.name))
        else:
            symbol = sym.Symbol(name)
            if not isinstance(orbitals, list):
                orbitals = [orbitals]

            # if symbol in self._species:
            #     s = self._species[symbol].create_sub_species(suffix=site.name)
            #     symbol = sym.Symbol(s.name)
            # else:
            s = Species(name, orbitals, site=site)
            # for name, new_species in s.sub_species.items():
            #     self._species[sym.Symbol(name)] = new_species

        self._species[symbol] = s
        return symbol

    def species_from_symbol(self, key: sym.Symbol) -> Species:
        return self._species[key]

    @property
    def sub_species_dict(self):
        """
        :return: a dictionary from all parent species to a list of its subspecies
        """
        d = {}
        for k, v in self._species.items():
            for n, sub_s in v.sub_species.items():
                if sym.Symbol(n) in self._species:
                    if k in d:
                        d[k].append(sym.Symbol(n))
                    else:
                        d[k] = [sym.Symbol(n)]
            if k in d and d[k]:
                d[k].append(k)
        return d


    def __str__(self):
        return str(self._species)  # TODO

    def __repr__(self):
        pass  # TODO

    sp = make_species
    __getitem__ = species_from_symbol
