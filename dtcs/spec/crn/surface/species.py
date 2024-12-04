"""TODO(Andrew)

"""

from __future__ import annotations
from typing import List, Union, Tuple

import itertools

import monty.json
import sympy as sym

import dtcs
from dtcs.common import util
from dtcs.spec.crn.surface.marker import Marker
from dtcs.spec.species import Species, SpeciesManager

# from dtcs.spec.xps import XPSSpecies, XPSSpeciesManager, XPSOrbital as Orbital


class SurfaceSpecies(Species):
    """A chemical species with a name and orbitals.

    Attributes:
        name: The name of the species.
        orbitals: A list of Orbitals.
    """

    def __init__(
        self,
        name: str,
        color: Union[Tuple[int], List[int], str] = None,
        parent=None,
        site: sym.Symbol = None,
        include_sub_species: bool = False,
        size: int = 1,
        is_gas: bool = False,
        is_visible: bool = False,
        **kwargs,
    ):
        super().__init__(
            name=name,
            color=color,
            **kwargs,
        )
        self.name = name
        self.is_gas = is_gas
        if not self.is_gas:
            if site is None:
                raise TypeError("Must supply site for non-gaseous species.")
            self.site = site.name
            self.size = size

        self.is_visible = is_visible

        # TODO(Andrew): Old stuff following
        # self.name_without_suffix = name

        # self.parent = parent   # The parent of the species
        # self.sub_species = {}   # The dictionary of sub species from name to the object

        # self.include_sub_species = include_sub_species
        # if include_sub_species and self.site and self.site != Site.default:

    @util.depreciate
    def create_sub_species(
        self,
        suffix: str = "",
        color: Union[Tuple[int], List[int], str] = "",
        entire_name: str = "",
        orbitals: Union[Orbital, List[Orbital]] = None,
        site: Site = None,
    ):
        if not entire_name:
            if site:
                suffix = site.name
            elif not suffix:
                suffix = f"sub_{len(self.sub_species)}"
            entire_name = f"{self.name_without_suffix}_{suffix}"

            # If a non-default site is added, add a new species
            self.name = f"{self.name_without_suffix}_{Site.default}"

        elif suffix:
            raise Exception(
                f"Both suffix={suffix} and entire_name {entire_name} provided to "
                + f"create_sub_species_function"
            )
        if not color:
            color = self.color
        if orbitals is None:
            orbitals = self.orbitals

        sub_species = SurfaceSpecies(
            entire_name,
            orbitals,
            color=color,
            parent=self,
            site=site,
            include_sub_species=False,
        )
        sub_species.name_without_suffix = self.name_without_suffix

        if entire_name in self.sub_species:
            raise Exception(
                f"species {self.name} already has sub species {repr(self.sub_species[entire_name])} "
                + f"with name {entire_name}."
            )
        else:
            self.sub_species[entire_name] = sub_species
        return sub_species

    @util.depreciate
    def sub_species_by_name(self, name: str):
        """
        :param name: name of the species
        :return: the species, if found
        """
        if name in self.sub_species:
            return self.sub_species[name]
        else:
            raise Exception(
                f"SurfaceSpecies {self.name} doesn't have a sub-species with name {name}."
            )

    @util.depreciate
    def sub_species_by_suffix(self, suffix: str):
        """
        :param suffix: the suffix of the species
        :return: the sub_species with the suffix, if found
        """
        return self.sub_species_by_name(self.sub_species_name(suffix))

    @util.depreciate
    def sub_species_name(self, suffix: str):
        """
        :param suffix: suffix of the sub_species
        :return: the sub_species's name
        """
        return f"{self.name_without_suffix}_{suffix}"

    # def __str__(self):
    #     orbitals = [str(orbital) for orbital in self.orbitals]
    #     if self.color:
    #         return f'{self.name}, orbitals: {orbitals}, color: {str(self.color)}, size: {str(self.size)}'
    #     return f'{self.name}, orbitals: {orbitals}, size: {str(self.size)}'
    #
    # def __repr__(self):
    #     return f'{self.__class__.__name__}(name={repr(self.name)}, ' \
    #            f'orbitals={repr(self.orbitals)}' + f'color={repr(self.color)}, size={self.size})'


class LegacySurfaceSpeciesManager(SpeciesManager):
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

    _species_cls = SurfaceSpecies

    @util.depreciate
    def __init__(self, *args, name="", **kwargs):
        name = name or f"Collection of {self._species_cls.__name__}"
        super().__init__(*args, name=name, **kwargs)

        self._species = {sym.Symbol(ele): ele for ele in self.elements}

        self.be_name_dict = {}  # What was this for ?
        self.default_surface_name = None
        self._markers = {}  # And what do these do?

    def add_species(self, species) -> sym.Symbol:
        s = super().add_species(species)
        self._species = {sym.Symbol(ele.name): ele for ele in self.elements}
        return s

    #
    # @property
    # def species(self) -> List[sym.Symbol]:
    #     """All of the species."""
    #     return list(self._species.keys())
    #
    # def make_species(self, name: Union[str, sym.Symbol],
    #                  orbitals: Union[Orbital, List[Orbital]] = None,
    #                  site: Site = None,
    #                  sub_species_name: str = "",
    #                  size=1,
    #                  color='', ) -> sym.Symbol:
    #     """
    #     Makes a sym.Symbol and a corresponding Species and keeps track of their correspondence.
    #     Keeps track of their correspondence.
    #     Orbitals can be either a list of orbitals or just one orbital
    #     TODO: use in the backend to consult for which products shall take two spots.
    #
    #     Args:
    #     name: The name of the new species and of the symbol.
    #     orbitals: The Orbitals of the species. Can be an Orbital or a list of Orbitals, just to be nice
    #     Returns:
    #     The sym.Symbol corresponding to the new Species.
    #     """
    #     if isinstance(name, sym.Symbol):
    #         parent = self.species_from_symbol(name)
    #
    #         if sub_species_name:
    #             name = sub_species_name
    #             if site:
    #                 name += site.name
    #             s = parent.create_sub_species(entire_name=name, site=site)
    #             symbol = sym.Symbol(name)
    #         else:
    #             s = parent.create_sub_species(site=site)
    #             symbol = sym.Symbol(parent.sub_species_name(site.name))
    #
    #             default_site_sym = sym.Symbol(parent.sub_species_name(Site.default))
    #             if default_site_sym not in self._species:
    #                 # TODO: change site.default to default name, as appropriate
    #                 self._species[default_site_sym] = parent.create_sub_species(site=Site(Site.default, parent.site))
    #     else:
    #         symbol = sym.Symbol(name)
    #         if orbitals is None:
    #             orbitals = []
    #         elif not isinstance(orbitals, list):
    #             orbitals = [orbitals]
    #
    #         # if symbol in self._species:
    #         #     s = self._species[symbol].create_sub_species(suffix=site.name)
    #         #     symbol = sym.Symbol(s.name)
    #         # else:
    #         # TODO: size for other occasions
    #         s = SurfaceSpecies(name, orbitals, site=site, size=size, color=color)
    #         # for name, new_species in s.sub_species.items():
    #         #     self._species[sym.Symbol(name)] = new_species
    #
    #     self._species[symbol] = s
    #     return symbol

    @util.depreciate
    def drop_marker(
        self, species_symbol: Union[sym.Symbol, Site], marker_name: str, color: str = ""
    ):
        """

        It's possible to mark two different species under the same name.
        :param species_symbol: symbol of the species being marked;
        :param name: name of the marker
        :return:
        """
        if marker_name in self._markers:
            for marker in self._markers[marker_name]:
                if marker.sm.name == str(species_symbol):
                    return marker
        new_marker = Marker(
            species_symbol.name,
            marker_name,
            species_symbol=sym.Symbol(species_symbol.name),
            color=color,
        )
        if marker_name in self._markers:
            self._markers[marker_name].append(new_marker)
        else:
            self._markers[marker_name] = [new_marker]
        return new_marker

    @util.depreciate
    def get_markers(self, name):
        return self._markers[name]

    @util.depreciate
    def get_marker_names(self):
        return self._markers.keys()

    @util.depreciate
    def get_all_markers(self):
        all_markers = []
        for markers in self._markers.values():
            all_markers.extend(markers)
        return all_markers

    @util.depreciate
    def name_be(
        self, name: str, value: float, orbital_name: str = "1S", color=""
    ) -> None:
        """
        name: the name for the binding energy
        value: numerical value of the binding energy
        """
        self.be_name_dict[value] = name
        # add a placeholder species
        self.make_species(name, [Orbital(orbital_name, value)], color=color)

    # def species_from_symbol(self, key: sym.Symbol) -> SurfaceSpecies:
    #     return self._species[key]

    def is_gas(self, species: Union[sym.Symbol, str]) -> bool:
        """
        Check if a species is gas.

        Rule: any species whose name does not end with "*" is a gas.

        :param species: a Sympy symbol referring to a gas species;
        :return: whether the symbol is a gas.
        """
        if isinstance(species, str):
            species = sym.Symbol(species)
        if self.symbol_in_sm(key=species) and self.species_from_symbol(species):
            return species.name.endswith("_g")
        return False

    def is_surface(self, species: Union[sym.Symbol, str]) -> bool:
        """
        Check if a species represents a surface.

        Rule: surface starts with "*".

        :param species: a Sympy symbol referring to a surface;
        :return: whether the symbol is a gas.
        """
        if isinstance(species, str):
            species = sym.Symbol(species)
        if self.symbol_in_sm(key=species) and self.species_from_symbol(species):
            return species.name.startswith("*")
        return False

    # def symbol_in_sm(self, key: sym.Symbol) -> bool:
    #     return key in self._species

    @property
    @util.depreciate
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
                #     all_default

            if k in d and d[k]:
                d[k].append(k)
        return d

    @property
    def to_sum_dict(self):
        """
        Provide users with dictionary to represent the summation relationships between species.

        :return: a list of names from species to sub-species and from binding energy to a list of species with the same binding energy.
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
                raise Exception(
                    "Our system currently doesn't support multiple orbitals"
                )
            be = species.orbitals[0].binding_energy
            if be in self.be_name_dict:
                be_name_symbol = sym.Symbol(self.be_name_dict[be])
                if be_name_symbol not in d:
                    d[be_name_symbol] = []
                d[be_name_symbol].append(sy)
                included_species_symbols.add(sy)
        return d

    # def symbol_from_name(self, name: str) -> sym.Symbol:
    #     """Gets the symbol for the given name, if it is a species.
    #
    #     Args:
    #         name: The name of the species you want to get the symbol for.
    #
    #     Raises:
    #         KeyError: If that's not a name for any species.
    #
    #     Returns:
    #         A sym.Symbol corresponding to a species in the species manager.
    #         sm.species_from_symbol(sm.symbol_from_name('name') will work
    #         unless you are supposed to get a KeyError.
    #     """
    #     symbol = sym.Symbol(name)
    #     if symbol in self._species:
    #         return symbol
    #     else:
    #         raise KeyError(f'Name {name} corresponds to no species.')
    #
    # def as_dict(self) -> dict:
    #     """Return a MSON-serializable dict representation."""
    #     d = {
    #         '@module': self.__class__.__module__,
    #         '@class': self.__class__.__name__,
    #         '@version': dtcs.__version__,  # TODO: Better way to do this?
    #         'species': {}
    #     }
    #
    #     for symbol, species in self._species.items():
    #         d['species'][symbol.name] = species.as_dict()
    #
    #     return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        decode = monty.json.MontyDecoder().process_decoded

        species_dict = {}
        for name, species in d["species"].items():
            species_dict[sym.Symbol(name)] = decode(species)
        d["species"] = species_dict

        return cls(**d)

    def __contains__(self, item: str):
        return item in self._species

    def __str__(self):
        s = self.__class__.__name__ + " with species:\n"
        for species in self._species.values():
            species_lines = str(species).splitlines()
            s += "".join([f"\t{line}\n" for line in species_lines])
        return s

    def __repr__(self):
        return f"{self.__class__.__name__}(species={repr(self._species)})"

    # sp = make_species
    # get = symbol_from_name
    # __getitem__ = species_from_symbol


class SurfaceSpeciesManager(SpeciesManager):

    _species_cls = SurfaceSpecies
