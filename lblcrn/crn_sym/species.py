"""
CRN - species.py:

Structures for the manipulation of chemical speices.

Credits:
Dr. Jin Qian, Domas Buracas, Ye Wang, Andrew Bogdan, Rithvik Panchapakesan
"""

from typing import List, Union

import monty.json
import sympy as sym


class Orbital(monty.json.MSONable):
    """
    An orbital in a species, this is essentially a named tuple, it's a class for readability purposes.
    """

    def __init__(self, name: str, binding_energy: float, splitting: float = 1):
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


class Species(monty.json.MSONable):
    """
    A chemical species with a name and orbitals, which are triples of (orbital name, binding energy, proportion)
    """

    def __init__(self, name: str, orbitals: List[Orbital]):
        self.name = name
        self.orbitals = orbitals

    def __str__(self):
        return self.name + ", orbitals: " + str(self.orbitals)

    def __repr__(self):
        return 'Species(name=' + self.name + ', orbitals=' + repr(self.orbitals) + ')'


class SpeciesManager(monty.json.MSONable):
    """
    A smart wrapper of a dictionary {sym.Symbol: Species} for the purpose of keeping track of
    which symbols correspond to which speices.

    You can create symbols/species pairs with SpeciesManager.sp and access them with SpeciesManager[],
    which forward to the more verbosely-named make_species and species_from_symbol
    """

    def __init__(self):
        self._species = {} # As of current, initializes empty

    def make_species(self, name: str, orbitals: Union[Orbital, List[Orbital]]) -> sym.Symbol:
        """
        Makes a sym.Symbol and a corresponding Species and keeps track of their correspondence.
        Returns the symbol.

        Orbitals can be either a list of orbitals or just one orbital, just to be kind.
        """
        symbol = sym.Symbol(name)

        if not isinstance(orbitals, list):
            orbitals = [orbitals]

        self._species[symbol] = Species(name, orbitals)
        return symbol

    def species_from_symbol(self, key: sym.Symbol) -> Species:
        return self._species[key]

    def symbol_from_name(self, name: str) -> sym.Symbol:
        """Gets the symbol for the given name, if it is a species.

        Args:
            name: The name of the species you want to get the symbol for.

        Raises:
            KeyError: If that's not a name for any species.

        Returns:
            A sym.Symbol corresponding to a species in the species manager.
            sm.species_from_symbol(sm.symbol_from_name('name') will work
            unless you are supposed to get a KeyError.
        """
        symbol = sym.Symbol(name)
        if symbol in self._species:
            return symbol
        else:
            raise KeyError(f'Name {name} corresponds to no species.')

    def __str__(self):
        return str(self._species)  # TODO

    def __repr__(self):
        pass  # TODO

    sp = make_species
    get = symbol_from_name
    __getitem__ = species_from_symbol
