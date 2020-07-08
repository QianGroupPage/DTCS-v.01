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

import monty.json
import sympy as sym
from lblcrn.common import color_to_RGB
from lblcrn.crn_sym.surface import Site
from typing import List, Tuple, Union

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

    def __init__(self, name: str, orbitals: List[Orbital], color: Union[Tuple[int], List[int], str] = None,
                 parent=None, site: Site = None, include_sub_species: bool=True, size: int = 1):
        self.name = name
        self.orbitals = orbitals
        if color:
            self.color = color_to_RGB(color)
        else:
            self.color = color
        self.parent = parent   # The parent of the species
        self.sub_species = {}   # The dictionary of sub species from name to the object

        self.site = site
        self.size = size
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
        orbitals = [str(orbital) for orbital in self.orbitals]
        if self.color:
            return f'{self.name}, orbitals: {orbitals}, color: {str(self.color)}, size: {str(self.size)}'
        return f'{self.name}, orbitals: {orbitals}, size: {str(self.size)}'

    def __repr__(self):
        return f'{self.__class__.__name__}(name={repr(self.name)}, ' \
               f'orbitals={repr(self.orbitals)}' + f'color={repr(self.color)}, size={self.size})'


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

    def make_species(self, name: Union[str, sym.Symbol], orbitals: Union[Orbital, List[Orbital]] = None,
                     site: Site = None, sub_species_name: str = "", size=1) -> sym.Symbol:
        """
        Makes a sym.Symbol and a corresponding Species and keeps track of their correspondence.

        Keeps track of their correspondence.
        Orbitals can be either a list of orbitals or just one orbital

        Args:
            name: The name of the new species and of the symbol.
            orbitals: The Orbitals of the species. Can be an Orbital or a list
                of Orbitals, just to be nice

        Returns:
             The sym.Symbol corresponding to the new Species.
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
            # TODO: size for other occasions
            s = Species(name, orbitals, site=site, size=size)
            # for name, new_species in s.sub_species.items():
            #     self._species[sym.Symbol(name)] = new_species

        self._species[symbol] = s
        return symbol

    def species_from_symbol(self, key: sym.Symbol) -> Species:
        return self._species[key]

    @property
    def all_species(self):
        return set(self._species.values())

    @property
    def large_species(self):
        res = set()
        for s in self.all_species:
            if s.size > 1:
                res.add(s)
        return res

    @property
    def all_symbols(self):
        return sorted(list(self._species.keys()), key=lambda s: str(s))

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
