"""TODO
"""

from lblcrn.spec import Spec
from lblcrn.spec.species import Species

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'


class XPSOrbital(Spec):
    """ TODO(Andrew): This is old

    An orbital in a species.

    This isn't actually a whole orbital. If you want to represent an orbital
    with splitting, you represent it with two orbitals, each with their own
    splitting coefficient.

    Attributes:
        name: The name of the orbital, e.g. 1s, 2p-1/2
        binding_energy: The binding energy of the orbital
        splitting: The splitting coefficient, these should sum to one.
        # TODO: Add the rest
    """

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

    def __init__(self, name: str, binding_energy: float, **kwargs):
        super().__init__(**kwargs)
        self.name = name
        self.binding_energy = binding_energy

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

    _schema = [
        'relax_vis',
        'xps_vis',
    ]

    _default = {
        'orbitals': list,
    }

    def __init__(self, name: str, **kwargs):
        """TODO"""
        if not ('orbitals' in kwargs or 'structure' in kwargs):
            raise TypeError('Either orbitals or structure required.')
        super().__init__(name=name, **kwargs)

    def __str__(self):
        # TODO
        orbitals = [str(orbital) for orbital in self.orbitals]
        if self.color:
            return f'{self.name}, orbitals: {orbitals}, color: {str(self.color)}, size: {str(self.size)}'
        return f'{self.name}, orbitals: {orbitals}, size: {str(self.size)}'

    def __repr__(self):
        # TODO
        return f'{self.__class__.__name__}(name={repr(self.name)}, ' \
               f'orbitals={repr(self.orbitals)}' + f'color={repr(self.color)}, size={self.size})'


# TODO: For backwards compatibility
Orbital = XPSOrbital
