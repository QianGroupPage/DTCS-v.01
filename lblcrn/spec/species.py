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

from __future__ import annotations
from typing import List, Union

import sympy as sym

from lblcrn.common import util
from lblcrn.common.colors import color_map
from lblcrn.spec.spec_abc import Spec, SpecCollection


class Species(Spec):
    """TODO"""

    def __init__(self, name, color='', **kwargs):
        """TODO"""
        super().__init__(name=name, **kwargs)
        if color:
            # Add color to the spec so it serializes
            self.spec['color'] = color
            # Use the global color tracking
            color_map[self.name] = color

    @property
    def color(self):
        return color_map[self.name]


class SpeciesManager(SpecCollection):
    """A smart wrapper of a dictionary {sym.Symbol: Species}. TODO(Andrew)

    Exists for the purpose of keeping track of which symbols correspond to
    which species.

    You can create symbols/species pairs with SpeciesManager.sp and access
    them with SpeciesManager[], which forward to the more verbosely-named
    make_species and species_from_symbol.

    If you need to, you can get a symbol which corresponds to a species with
    SpeciesManager.get, which forwards to symbol_from_name. This is useful if,
    for example you loaded the SpeciesManager from a file.
    """

    _species_cls = Species

    def __init__(self, *species, name=None, **kwargs):
        name = name or f'Collection of {self._species_cls.__name__}'
        super().__init__(*species, name=name, **kwargs)

    # --- Helpful Accessors ---------------------------------------------------
    @property
    def names(self) -> List[str]:  # TODO: Move or remove?
        return [species.name for species in self.elements]

    @property
    def species(self) -> List[Species]:  # TODO: Move or remove?
        """All of the species."""
        return list(self.elements)

    @property
    def symbols(self) -> List[sym.Symbol]:  # TODO: Move or remove?
        """Symbols for all of the species. Fixed order."""
        return [sym.Symbol(name) for name in self.names]

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
        if name in self._names:
            return sym.Symbol(name)
        else:
            raise KeyError(f'Name {name} corresponds to no species.')

    def add_species(self, species: Species) -> sym.Symbol:
        """TODO"""
        if not species.name:
            raise ValueError('The species is missing a name.')
        if species.name in self._names:
            raise ValueError(f'Species {species.name} is already in this'
                             f'SpeciesCollection.')

        self.append(species)
        return sym.Symbol(species.name)

    def make_species(self, *args, **kwargs) -> sym.Symbol:
        """Makes a sym.Symbol and a corresponding Species.
        Keeps track of their correspondence.

        TODO(Andrew): Args
        Returns:
            The sym.Symbol corresponding to the new Species.
        """
        return self.add_species(self._species_cls(*args, **kwargs))

    # --- Override Magic to accept sym.Symbols too ----------------------------
    def __getitem__(self, key: Union[str, sym.Symbol]) -> Species:
        return super().__getitem__(util.symbol_to_name(key))

    def __setitem__(self, key: Union[str, sym.Symbol], value: Species):
        return super().__setitem__(util.symbol_to_name(key), value)

    def __delitem__(self, key: Union[str, sym.Symbol]):
        return super().__delitem__(util.symbol_to_name(key))

    def __contains__(self, item: Union[str, sym.Symbol]) -> bool:
        """Convert symbols to strings and then """
        return super().__contains__(util.symbol_to_name(item))

    # --- Other Magic ---------------------------------------------------------
    def __str__(self):
        s = self.__class__.__name__ + ' with species:\n'
        for species in self.elements:
            species_lines = str(species).splitlines()
            s += ''.join([f'\t{line}\n' for line in species_lines])
        return s

    def __repr__(self):
        return f'{self.__class__.__name__}(species={repr(self.elements)})'

    # Useful shorthand
    add = add_species
    get = symbol_from_name
    sp = make_species

    @util.depreciate
    def species_from_symbol(self, key: sym.Symbol) -> Species:  # TODO: Depreciate
        return self[key]

    @util.depreciate
    def symbol_in_sm(self, key: sym.Symbol) -> bool:  # TODO: Depreciate
        return key in self._names

    @property
    @util.depreciate
    def symbols_ordering(self):  # TODO: Depreciate
        # TODO: This used to be sorted(self._species, key=lambda s: str(s)),
        #  check if the alphabetical order was actually used.
        return list(self.symbols)

    @property
    @util.depreciate
    def all_species(self):  # TODO: Depreciate
        return set(self.elements)

    @property
    @util.depreciate
    def all_symbols(self):  # TODO: Depreciate
        return self.symbols