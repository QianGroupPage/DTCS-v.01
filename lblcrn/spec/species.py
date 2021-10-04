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

import itertools
from typing import List, Union

import sympy as sym

from lblcrn.spec import Spec
from lblcrn.common import util


_COLORS = itertools.cycle(['red', 'green', 'orange', 'blue', 'purple', 'pink',
                           'yellow', 'gray', 'cyan'])


class Species(Spec):
    """TODO"""

    _schema = [
        'name',
        'structure',
    ]

    _default = {
        'color': '',
    }

    def __init__(self, name, **kwargs):
        """TODO"""
        super().__init__(**kwargs)
        self.name: str = name


class SpeciesCollection(Spec):
    """TODO"""

    @property
    def names(self) -> List[str]:
        return list(self.spec.keys())

    @property
    def species(self) -> List[Species]:
        """All of the species."""
        return list(self.spec.values())

    @property
    def symbols(self) -> List[sym.Symbol]:
        """Symbols for all of the species."""
        return [sym.Symbol(name) for name in self.spec]

    def add_species(self, species: Species) -> sym.Symbol:
        """TODO"""
        if not species.name:
            raise ValueError('The species is missing a name.')
        if species.name in self.spec:
            raise ValueError(f'Species {species.name} is already in this'
                             f'SpeciesCollection.')

        self.spec[species.name] = species
        return sym.Symbol(species.name)

    def make_species(self, **kwargs) -> sym.Symbol:
        """TODO"""
        return self.add_species(Species(**kwargs))

    # Mapper magic methods take both strings and sym.Symbols
    def __getitem__(self, key: Union[str, sym.Symbol]) -> Species:
        return super().__getitem__(util.symbol_to_name(key))

    def __setitem__(self, key: Union[str, sym.Symbol], value: Species):
        # Ensure that keys and species names match
        if key != value.name:
            raise ValueError('Keys and species names must match.')
        return super().__setitem__(util.symbol_to_name(key), value)

    def __delitem__(self, key: Union[str, sym.Symbol]):
        return super().__delitem__(util.symbol_to_name(key))

    def __contains__(self, item: Union[str, sym.Symbol]) -> bool:
        """Convert symbols to strings and then """
        return super().__contains__(util.symbol_to_name(item))

    def __getattr__(self, name):
        """Don't allow accessing species like attributes."""
        raise AttributeError(f'\'{self.__class__.__name__}\' object'
                             f'has no attribute \'{name}\'')

    # Don't allow accessing species like attributes
    __setattr__ = object.__setattr__
    __delattr__ = object.__delattr__

    # Useful shorthand
    add = add_species
    sp = make_species
