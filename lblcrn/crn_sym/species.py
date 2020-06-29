"""Classes for defining chemical species.

Ideally we'd transition to using pymatgen's stuff, but before then, this is
good.

Exports:
    Orbital: An orbital in a species.
    Species: A chemical species with a name and orbitals.
    SpeciesManager: A smart wrapper of a dictionary {sym.Symbol: Species}.

Usage:
    # Create and print a new species X.
    sm = SpeciesManager()
    x = sm.sp('X', Orbital('1s', 340))
    print(sm[x])
"""

from typing import List, Union

import lblcrn
import monty.json
import sympy as sym


class Orbital(monty.json.MSONable):
    """An orbital in a species.

    This isn't actually a whole orbital. If you want to represent an orbital
    with splitting, you represent it with two orbitals, each with their own
    splitting coefficient.

    Attributes:
        name: The name of the orbital, e.g. 1s, 2p-1/2
        binding_energy: The binding energy of the orbital
        splitting: The splitting coefficient, these should sum to one.
    """

    def __init__(self, name: str, binding_energy: float, splitting: float = 1):
        self.name = name
        self.binding_energy = binding_energy
        self.splitting = splitting

    def __str__(self):
        if self.splitting == 1:
            return f'{self.name} @ {self.binding_energy} eV'
        else:
            return f'{self.name} @ {self.binding_energy} eV, ' \
                   f'splitting {self.splitting}'

    def __repr__(self):
        if self.splitting == 1:
            return f'{self.__class__.__name__}(name={repr(self.name)}, ' \
                   f'binding_energy={repr(self.binding_energy)})'
        else:
            return f'{self.__class__.__name__}(name={repr(self.name)}, ' \
                   f'binding_energy={repr(self.binding_energy)}, ' \
                   f'splitting={repr(self.splitting)})'


class Species(monty.json.MSONable):
    """A chemical species with a name and orbitals.

    Attributes:
        name: The name of the species.
        orbitals: A list of Orbitals.
    """

    def __init__(self, name: str, orbitals: List[Orbital]):
        self.name = name
        self.orbitals = orbitals

    def __str__(self):
        orbitals = [str(orbital) for orbital in self.orbitals]
        return f'{self.name}, orbitals: {orbitals}'

    def __repr__(self):
        return f'{self.__class__.__name__}(name={repr(self.name)}, ' \
               f'orbitals={repr(self.orbitals)})'


class SpeciesManager(monty.json.MSONable):
    """A smart wrapper of a dictionary {sym.Symbol: Species}.

    Exists for the purpose of keeping track of which symbols correspond to
    which speices.

    You can create symbols/species pairs with SpeciesManager.sp and access
    them with SpeciesManager[], which forward to the more verbosely-named
    make_species and species_from_symbol.

    If you need to, you can get a symbol which corresponds to a species with
    SpeciesManager.get, which forwards to symbol_from_name. This is useful if,
    for example you loaded the SpeciesManager from a file.
    """

    def __init__(self, species: dict = None):
        if species:
            self._species = species
        else:
            self._species = {}

    def make_species(self, name: str,
                     orbitals: Union[Orbital, List[Orbital]]) -> sym.Symbol:
        """Makes a sym.Symbol and a corresponding Species

        Keeps track of their correspondence.
        Orbitals can be either a list of orbitals or just one orbital

        Args:
            name: The name of the new species and of the symbol.
            orbitals: The Orbitals of the species. Can be an Orbital or a list
                of Orbitals, just to be nice

        Returns:
             The sym.Symbol corresponding to the new Species.
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

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = {
            '@module': self.__class__.__module__,
            '@class': self.__class__.__name__,
            '@version': lblcrn.__version__,  # TODO: Better way to do this?
            'species': {}
        }

        for symbol, species in self._species.items():
            d['species'][symbol.name] = species.as_dict()

        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        decode = monty.json.MontyDecoder().process_decoded

        species_dict = {}
        for name, species in d['species'].items():
            species_dict[sym.Symbol(name)] = decode(species)
        d['species'] = species_dict

        return cls(**d)

    def __str__(self):
        s = self.__class__.__name__ + ' with species:\n'
        for species in self._species.values():
            species_lines = str(species).splitlines()
            s += ''.join([f'\t{line}\n' for line in species_lines])
        return s

    def __repr__(self):
        return f'{self.__class__.__name__}(species={repr(self._species)})'

    sp = make_species
    get = symbol_from_name
    __getitem__ = species_from_symbol
