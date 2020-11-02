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
from typing import List, Optional, Tuple, Union

import monty.json
import sympy as sym

import lblcrn
from lblcrn.crn_sym.surface import Site

_COLORS = itertools.cycle(['red', 'green', 'orange', 'blue', 'purple', 'pink',
                           'yellow', 'gray', 'cyan'])


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
        self.name_without_suffix = name

        self.orbitals = orbitals
        self.color = color or next(_COLORS)
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
            entire_name = f"{self.name_without_suffix}_{suffix}"

            # If a non-default site is added, add a new species
            self.name = f"{self.name_without_suffix}_{Site.default}"

        elif suffix:
            raise Exception(f"Both suffix={suffix} and entire_name {entire_name} provided to " +
                            f"create_sub_species_function")
        if not color:
            color = self.color
        if orbitals is None:
            orbitals = self.orbitals

        sub_species = Species(entire_name, orbitals, color=color, parent=self, site=site, include_sub_species=False)
        sub_species.name_without_suffix = self.name_without_suffix

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
        return f"{self.name_without_suffix}_{suffix}"

    def __str__(self):
        orbitals = [str(orbital) for orbital in self.orbitals]
        if self.color:
            return f'{self.name}, orbitals: {orbitals}, color: {str(self.color)}, size: {str(self.size)}'
        return f'{self.name}, orbitals: {orbitals}, size: {str(self.size)}'

    def __repr__(self):
        return f'{self.__class__.__name__}(name={repr(self.name)}, ' \
               f'orbitals={repr(self.orbitals)}' + f'color={repr(self.color)}, size={self.size})'


class Marker(monty.json.MSONable):
    """
    A marker object to mark occurances of a given species with a different name.
    """
    def __init__(self, species: str, name: str, species_symbol: sym.Symbol = None, color: str = ""):
        self.species = species

        self.species_symbol = species_symbol if species_symbol else sym.Symbol(species.name)
        self.name = name
        self.initial_count = 0
        self.color = color

    def __str__(self):
        return f'Marker with name {self.name} for species {repr(self.species)}'

    def __repr__(self):
        return f'{self.__class__.__name__}(name={repr(self.name)}, species={repr(self.species)}'


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

    def __init__(self, species: dict = None, colors: dict = None):
        if species:
            self._species = species
        else:
            self._species = {}

        self.be_name_dict = {}
        self.default_surface_name = None

        self._markers = {}

    @property
    def species(self) -> List[sym.Symbol]:
        """All of the species."""
        return list(self._species.keys())

    def make_species(self, name: Union[str, sym.Symbol],
                     orbitals: Union[Orbital, List[Orbital]] = None,
                     site: Site = None,
                     sub_species_name: str = "",
                     size=1,
                     color='', ) -> sym.Symbol:
        """
        Makes a sym.Symbol and a corresponding Species and keeps track of their correspondence.
        Keeps track of their correspondence.
        Orbitals can be either a list of orbitals or just one orbital
        TODO: use in the backend to consult for which products shall take two spots.
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

                default_site_sym = sym.Symbol(parent.sub_species_name(Site.default))
                if default_site_sym not in self._species:
                    # TODO: change site.default to default name, as appropriate
                    self._species[default_site_sym] = parent.create_sub_species(site=Site(Site.default, parent.site))
        else:
            symbol = sym.Symbol(name)
            if orbitals is None:
                orbitals = []
            elif not isinstance(orbitals, list):
                orbitals = [orbitals]

            # if symbol in self._species:
            #     s = self._species[symbol].create_sub_species(suffix=site.name)
            #     symbol = sym.Symbol(s.name)
            # else:
            # TODO: size for other occasions
            s = Species(name, orbitals, site=site, size=size, color=color)
            # for name, new_species in s.sub_species.items():
            #     self._species[sym.Symbol(name)] = new_species

        self._species[symbol] = s
        return symbol

    def drop_marker(self, species_symbol: Union[sym.Symbol, Site], marker_name: str, color: str=""):
        """

        It's possible to mark two different species under the same name.
        :param species_symbol: symbol of the species being marked;
        :param name: name of the marker
        :return:
        """
        if marker_name in self._markers:
            for marker in self._markers[marker_name]:
                if marker.species.name == str(species_symbol):
                    return marker
        new_marker = Marker(species_symbol.name, marker_name, species_symbol=sym.Symbol(species_symbol.name),
                            color=color)
        if marker_name in self._markers:
            self._markers[marker_name].append(new_marker)
        else:
            self._markers[marker_name] = [new_marker]
        return new_marker

    def get_markers(self, name):
        return self._markers[name]

    def get_marker_names(self):
        return self._markers.keys()

    def get_all_markers(self):
        all_markers = []
        for markers in self._markers.values():
            all_markers.extend(markers)
        return all_markers

    def name_be(self, name: str, value: float, orbital_name: str = "1S", color="") -> None:
        """
        name: the name for the binding energy
        value: numerical value of the binding energy
        """
        self.be_name_dict[value] = name
        # add a placeholder species
        self.make_species(name, [Orbital(orbital_name, value)], color=color)

    def species_from_symbol(self, key: sym.Symbol) -> Species:
        return self._species[key]

    def symbol_in_sm(self, key: sym.Symbol) -> bool:
        return key in self._species

    @property
    def symbols_ordering(self):
        return sorted(self._species, key=lambda s: str(s))

    def get_site_name(self, species_name: str) -> str:
        """
        Return the site name of the site where the species is supposed to be.
        """
        if sym.Symbol(species_name) in self._species:
            site = self.species_from_symbol(sym.Symbol(species_name)).site
            if site is None:
                return self.default_surface_name
            elif site.name == Site.default:
                return site.surface.name
            else:
                return site.name
        # the species_name refers to a site
        else:
            return species_name

    # TODO: test this new function
    @property
    def names_by_site_name(self):
        result = {}
        for species in self.all_species:
            site_name = self.get_site_name(species.name)
            if site_name not in result:
                result[site_name] = []
            result[site_name].append(species.name)
        return result

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
    def large_species_dict(self):
        """
        Return a dictionary. This is primarily for use inside the surface crn back-end.
        """
        return {species.name: species.size for species in self.large_species}

    @property
    def all_symbols(self):
        return sorted(list(self._species.keys()), key=lambda s: str(s))

    @property
    def aggregated_symbols(self):
        """
        1. For species with various sites, only the symbol containing the default name would be included.

        2. AFTER RULE 1 IS EXECUTED, symbols with the same binding energy shall be grouped together.
        """
        pass

    @property
    def backend_symbols(self):
        """
        Symbols corresponding to species used in the simulation backend.
        """
        pass

    @property
    def sub_species_dict(self):
        """
        :return: a dictionary from all parent species symbol to a list of its subspecies symbols
        """
        d = {}
        for k, v in self._species.items():
            all_default = False
            for n, sub_s in v.sub_species.items():
                if sym.Symbol(n) in self._species:
                    if k in d:
                        d[k].append(sym.Symbol(n))
                    else:
                        d[k] = [sym.Symbol(n)]

                # if sub_s.site.name == Site.default:
                #     all_defauly
                #

            if k in d and d[k]:
                d[k].append(k)
        return d

    @property
    def to_sum_dict(self):
        """
        Provide users with dictionary to represent the summation relationships between species.

        :return: a list of names from species to sub-species and from binding energy to a list of species
        with the same binding energy.
        """
        d = {}
        included_species_symbols = set()
        for species_symbol, l in self.sub_species_dict.items():
            d[species_symbol] = l
            for sub_s_symbol in l:
                included_species_symbols.add(sub_s_symbol)

        for sy in self.all_symbols:
            if sy in included_species_symbols:
                continue
            species = self.species_from_symbol(sy)
            if len(species.orbitals) != 1:
                # TODO
                raise Exception("Our system currently doesn't support multiple orbitals")
            be = species.orbitals[0].binding_energy
            if be in self.be_name_dict:
                be_name_symbol = sym.Symbol(self.be_name_dict[be])
                if be_name_symbol not in d:
                    d[be_name_symbol] = []
                d[be_name_symbol].append(sy)
                included_species_symbols.add(sy)
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
