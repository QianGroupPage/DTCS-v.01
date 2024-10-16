from typing import List
from dtcs.spec.species import Species


__author__ = "Rahul Khorana"


class IRSpecies(Species):
    """
    Address each species as a dictionary.
    Integrated into -> twin -> xps.py

    """

    def __init__(self, speciesList: List):
        self.sl = speciesList  # species list -> for each species is a dictionary
        # dictionary: (number -> location, intensity)
        # for each species -> generate standard spectra
        # units: IR -> transmittance; wavenumber cm^-1

    def plot_spectra(self):
        """plot the IR spectra -> list of separate graphs"""
        return
