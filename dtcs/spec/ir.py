from typing import List, Tuple

from dtcs.spec.species import Species, Spec
from dtcs.spec.xps import XPSSpeciesManager

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

__author__ = "Rahul Khorana"


class IRSignature(Spec):
    def __init__(self, name: str, peaks: dict, **kwargs):
        """
        TODO: Example of what is acceptable
        """
        super().__init__(name=name, **kwargs)
        try:
            assert isinstance(name, str)
            assert isinstance(peaks, dict) and len(peaks.keys()) == 2
        except:
            error = "Invalid input: name needs to be a string and peaks needs to be a dictionary \n"
            error += "peaks must contain two keys {'location','intensity'}"
            raise Exception(error)
        allowed_keys = set(["location", "intensity"])
        for k in peaks.keys():
            assert k in allowed_keys
        for v in peaks.values():
            try:
                assert isinstance(v, list) and len(v) > 0
                types = set([float, int, np.float32, np.float16, np.float64])
                for item in v:
                    assert type(item) in types
            except:
                error = "The lists must be nonempty and contain either: float, int, or np.float"
                raise Exception(error)
        for intensity in peaks["intensity"]:
            try:
                assert 0 <= intensity <= 1
            except:
                raise Exception("intensity must be between 0 and 1")
        self.peaks = peaks
        self.name = name

    def __str__(self):
        return f"name: {self.name} peaks: {self.peaks}"


class IRSpecies(Species):
    def __init__(
        self,
        name: str,
        ir_signature: Tuple[List[IRSignature] | IRSignature],
        is_visible: bool = False,
    ):
        """
        # species list -> for each species is a dictionary
        # dictionary: (number -> location, intensity)
        # for each species -> generate standard spectra
        # units: IR -> transmittance; wavenumber cm^-1
        TODO: Example of what is acceptable
        """
        assert isinstance(name, str)
        assert (isinstance(ir_signature, list) and len(ir_signature) > 0) or isinstance(
            ir_signature, IRSignature
        )
        assert isinstance(is_visible, bool)
        if isinstance(ir_signature, list):
            for i in range(len(ir_signature)):
                ith_species = ir_signature[i]
                assert isinstance(ith_species, IRSignature)
        if isinstance(ir_signature, IRSignature):
            ir_signature = [ir_signature]
        self.species_list = ir_signature
        self.is_visible = is_visible

    def plot_spectra(self):
        """plot the IR spectra -> list of separate graphs"""
        for i, species in enumerate(self.species_list):
            assert isinstance(species, IRSignature)
            path = f"spectra_plots/plot_{i}"
            self.plot_single_spectra(species.peaks, path)
        return

    def plot_single_spectra(self, spectra_data, path):
        wavenumbers = spectra_data["location"]
        transmittance = spectra_data["intensity"]
        sns.set_theme(style="whitegrid")
        sns.set_context("talk")
        plt.figure(figsize=(10, 6))
        sns.lineplot(x=wavenumbers, y=transmittance, linewidth=2.5)
        plt.gca().invert_xaxis()
        plt.xlabel("Wavenumber (cm$^{-1}$)", fontsize=14)
        plt.ylabel("Transmittance (%)", fontsize=14)
        plt.title(
            "Spectra Plot: Transmittance vs Wavenumber", fontsize=16, fontweight="bold"
        )
        plt.grid(True, which="both", axis="both", linestyle="--", linewidth=0.7)
        plt.tight_layout()
        plt.savefig(path)


if __name__ == "__main__":
    spectra_data_1 = {
        "location": [4000, 3500, 3000, 2500, 2000, 1500, 1000],
        "intensity": [95, 90, 85, 80, 75, 70, 65],
    }
    spectra_data_2 = {
        "location": [
            4000,
            3950,
            3900,
            3850,
            3800,
            3750,
            3700,
            3650,
            3600,
            3550,
            3500,
            3450,
            3400,
            3350,
            3300,
            3250,
            3200,
            3150,
            3100,
            3050,
            3000,
            2950,
            2900,
            2850,
            2800,
            2750,
            2700,
            2650,
            2600,
            2550,
            2500,
            2450,
            2400,
            2350,
            2300,
            2250,
            2200,
            2150,
            2100,
            2050,
            2000,
            1950,
            1900,
            1850,
            1800,
            1750,
            1700,
            1650,
            1600,
            1550,
            1500,
            1450,
            1400,
            1350,
            1300,
            1250,
            1200,
            1150,
            1100,
            1050,
            1000,
        ],
        "intensity": [
            95,
            94,
            93,
            91,
            88,
            86,
            85,
            83,
            80,
            78,
            76,
            75,
            72,
            70,
            69,
            66,
            64,
            63,
            61,
            60,
            59,
            57,
            55,
            53,
            52,
            51,
            50,
            49,
            48,
            47,
            45,
            43,
            42,
            40,
            39,
            37,
            36,
            34,
            32,
            31,
            30,
            29,
            28,
            26,
            25,
            23,
            21,
            20,
            19,
            18,
            16,
            15,
            14,
            13,
            12,
            11,
            10,
            9,
            8,
            7,
            6,
        ],
    }
    irs = IRSpecies(speciesList=[spectra_data_1, spectra_data_2])
    irs.plot_spectra()
