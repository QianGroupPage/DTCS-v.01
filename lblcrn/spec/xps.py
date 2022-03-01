"""TODO
"""

from __future__ import annotations
from typing import Optional, List, Union

import monty.json
from pymatgen.core.structure import Structure

from lblcrn.spec.spec_abc import Spec
from lblcrn.spec.species import Species, SpeciesManager
from lblcrn.spec.crn.surface.species import SurfaceSpecies

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'

class XPSOrbital(Spec):
    """An orbital in a species.

    This isn't actually a whole orbital. If you want to represent an orbital
    with splitting, you represent it with two orbitals, each with their own
    splitting coefficient.

    Attributes:
        name: The name of the orbital, e.g. 1s, 2p-1/2
        binding_energy: The binding energy of the orbital
        splitting: The splitting coefficient, these should sum to one.

    TODO(Andrew) things to add to this
    _schema = [
        'name',
        'element',
        'orbital',
        'binding_energy',
        'site_num',
    ]

    _default = {
        'splitting': 1,
        'is_surface': False,
    }
    """

    def __init__(self, name: str, binding_energy: float, splitting: float = 1, **kwargs):
        super().__init__(name=name, **kwargs)
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


class XPSSpecies(Species):
    """TODO"""

    def __init__(self,
                 name: str,
                 orbs_or_struct: Union[List, Structure] = None,
                 *,
                 orbitals: List[XPSOrbital] = None,
                 structure: Structure = None,
                 relax_vis=None,
                 xps_vis=None,
                 **kwargs):
        """TODO: Write docstring, check types, make it better, etc.
        """
        if orbitals:
            if isinstance(orbitals, XPSOrbital):
                orbitals = [orbitals]
            self.orbitals = orbitals
        elif isinstance(orbs_or_struct, (list, tuple)):
            self.orbitals = orbs_or_struct
        elif isinstance(orbs_or_struct, XPSOrbital):
            self.orbitals = [orbs_or_struct]
        elif isinstance(orbs_or_struct, (int, float)):
            self.orbitals = [XPSOrbital(name,
                                        binding_energy=orbs_or_struct)]

        if structure:
            self.structure = structure
        elif isinstance(orbs_or_struct, Structure):
            self.structure = orbs_or_struct

        if not (hasattr(self, 'orbitals') or hasattr(self, 'structure')):
            raise TypeError('Either orbitals or structure required.')

        if relax_vis: self.relax_vis = relax_vis
        if xps_vis: self.xps_vis = xps_vis
        Species.__init__(self, name=name, **kwargs)  # TODO(Andrew): Bodge

    def __str__(self):
        # TODO: will break if there's a structure, not orbitals
        orbitals = [str(orbital) for orbital in self.orbitals]
        if self.color:
            return f'{self.name}, orbitals: {orbitals}, color: {str(self.color)}'
        return f'{self.name}, orbitals: {orbitals}'

    def __repr__(self):
        # TODO: will break if there's a structure, not orbitals
        return f'{self.__class__.__name__}(name={repr(self.name)}, ' \
               f'orbitals={repr(self.orbitals)}' + f'color={repr(self.color)}'


class XPSSurfaceSpecies(XPSSpecies, SurfaceSpecies):
    def __init__(self,
                 name: str,
                 orbs_or_struct: Union[List, Structure] = None,
                 *,
                 orbitals: List[XPSOrbital] = None,
                 structure: Structure = None,
                 relax_vis=None,
                 xps_vis=None,

                 site=None,
                 size=None,
                 is_gas=False,
                 **kwargs):

        SurfaceSpecies.__init__(
            self,
            name=name,
            site=site,
            size=size,
            is_gas=is_gas,
        )

        XPSSpecies.__init__(
            self,
            name=name,
            orbs_or_struct=orbs_or_struct,
            orbitals=orbitals,
            structure=structure,
            relax_vis=relax_vis,
            xps_vis=xps_vis,
            **kwargs
        )


class XPSSpeciesManager(SpeciesManager):
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

    _species_cls = XPSSpecies

    # def name_be(self, name: str, value: float, orbital_name: str = "1S", color="") -> None:
    #     """
    #     name: the name for the binding energy
    #     value: numerical value of the binding energy
    #     """
    #     self.be_name_dict[value] = name
    #     # add a placeholder species
    #     self.make_species(name, [Orbital(orbital_name, value)], color=color)

    # TODO: These two functions are a really good idea, I removed them
    #  to clean up, but once everything has stabilized, I should move them from
    #  surface/species.py
    # def is_gas(self, species: Union[sym.Symbol, str]) -> bool:
    # def is_surface(self, species: Union[sym.Symbol, str]) -> bool:


# TODO: For backwards compatibility
Orbital = XPSOrbital
