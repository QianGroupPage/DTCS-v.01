"""
...
"""

# *** Libraries ***
import pkg_resources

import pandas as pd
import sympy as sym
from typing import Iterator, Tuple, Union

from .species import *

# *** Resources ***
binding_energy_csv = pkg_resources.resource_filename('lblcrn', 'resources/binding_energy_data.csv')

# *** Constants ***
BINDING_ENERGY_DATA = pd.read_csv(binding_energy_csv, index_col=0)
ELEMENTS = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 
            'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 
            'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 
            'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 
            'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 
            'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 
            'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 
            'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U']
ELEMENT_INDEX = list(range(1, 93))
ORBITALS = {'1s': ['1s'],
     '2s': ['2s'],
     '2p': ['2p1-2', '2p3-2'],
     '3s': ['3s'],
     '3p': ['3p1-2', '3p3-2'],
     '3d': ['3d3-2', '3d5-2'],
     '4s': ['4s'],
     '4p': ['4p1-2', '4p3-2'],
     '4d': ['4d3-2', '4d5-2'],
     '4f': ['4f5-2', '4f7-2'],
     '5s': ['5s'],
     '5p': ['5p1-2', '5p3-2'],
     '5d': ['5d3-2', '5d5-2'],
     '6s': ['6s'],
     '6p': ['6p1-2', '6p3-2']}
ORBITAL_SPLITS = {'1s': '1s',
     '2s': '2s',
     '2p1-2': '2p',
     '2p3-2': '2p',
     '3s': '3s',
     '3p1-2': '3p',
     '3p3-2': '3p',
     '3d3-2': '3d',
     '3d5-2': '3d',
     '4s': '4s',
     '4p1-2': '4p',
     '4p3-2': '4p',
     '4d3-2': '4d',
     '4d5-2': '4d',
     '4f5-2': '4f',
     '4f7-2': '4f',
     '5s': '5s',
     '5p1-2': '5p',
     '5p3-2': '5p',
     '5d3-2': '5d',
     '5d5-2': '5d',
     '6s': '6s',
     '6p1-2': '6p',
     '6p3-2': '6p'}

# *** Classes ***
#@pd.api.extensions.register_dataframe_accessor("be")
class BindingEnergyAccessor:
    """
    ...
    """
    
    def __init__(self, pandas_obj, sm):
        self._validate(pandas_obj)
        self._obj = pandas_obj
        self._sm = sm
    
    @staticmethod
    def _validate(obj):
        """FIXME"""
        pass
    
    def get_species(self, element, orbital_split):
        """
        Gets a Species of element at the orbital of orbital_split.
        """
        # Ensure that orbital is an orbtial, not just a split.
        if orbital_split in ORBITAL_SPLITS:
            orbital = ORBITAL_SPLITS[orbital_split]
        else:
            orbital = orbital_split
            
        # Collect the Orbitals
        orbitals = []
        for split in ORBITALS[orbital]:
            orbitals.append(Orbital(split, self._obj[split][element], 1)) #FIXME: Add splitting
            
        # Return the Species
        return self._sm.sp(element, orbitals, {}) #FIXME: Empty schedule
        
    def search(self, binding_energy: float, percent_margin: float = 0.01) -> Iterator[Tuple[str, str, int]]:
        """
        Yield all orbitals within percent_margin of binding_energy.
        """
        lower = (1 - percent_margin) * binding_energy
        upper = (1 + percent_margin) * binding_energy
        
        for element, orbital_split, binde in self:
            if binde >= lower and binde <= upper:
                yield element, orbital_split, binde
    
    def select(self, binding_energy: float, percent_margin: float = 0.01) -> sym.Symbol:
        """
        Go through orbitals within percent_margin of binding_energy, ask 
        with input() if that is the orbital you're looking for, and then 
        return the species's symbol, if you say yes, adding that species 
        to your SpeciesManager.
        
        Meant for use in the interpreter or in a notebook.
        """
        for element, orbital_split, binde in self.search(binding_energy, percent_margin):
            correct = input(f'Is {element} {orbital_split} at {binde} eV correct? (y/n): ')
            if 'y' in correct or 'Y' in correct:
                return self.get_species(element, orbital_split)
            
    def __getitem__(self, key: Union[int, str]) -> pd.DataFrame:
        """
        Gives rows by element, by orbital, or by orbital split.
        """
        if key in ORBITAL_SPLITS.keys():
            # Return the column of that orbital split
            return self._obj[[key,]]
        elif key in ORBITALS.keys():
            # Return the columns of all that orbital's splits
            cols = [col for col in bed.columns if col in ORBITALS[key]]
            return self._obj[cols]
        elif key in ELEMENTS:
            # Return the row of that element (by name)
            return self._obj.loc[[key,]]
        elif isinstance(key, int) and key < len(ELEMENTS):
            # Return the row of that element 
            return self._obj.loc[[ELEMENTS[key],]]
        raise KeyError(key)
        
    def __iter__(self) -> Iterator[Tuple[str, str, int]]:
        """
        Iterates over every non-NaN datapoint, yielding tuples of the 
        form (element, orbital split, binding energy).
        
        I suggest you use it like
        `for element, orbital_split, binding_energy in <BindingEnergyAccessor>:
            ...`.
        """
        for element in ELEMENTS:
            for orbital_split in ORBITAL_SPLITS.keys():
                if not pd.isna(self._obj[orbital_split][element]):
                    yield element, orbital_split, self._obj[orbital_split][element]
        
be = BindingEnergyAccessor(BINDING_ENERGY_DATA, SpeciesManager())