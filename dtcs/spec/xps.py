"""TODO
"""

from __future__ import annotations

import logging
import warnings
from typing import Optional, List, Union, Tuple, Dict, TypeVar

_logger = logging.getLogger(__name__)

import monty.json
import numpy as np
import pandas as pd
import sympy as sym

try:
    from pymatgen.core.structure import Structure
except ModuleNotFoundError:
    Structure = TypeVar("Structure")
    _logger.info("Didn't load module pymatgen")

from dtcs import _logger
from dtcs.sim.xps import simulate_xps

from dtcs.spec.spec_abc import Spec
from dtcs.spec.species import Species, SpeciesManager
from dtcs.spec.crn.surface.species import SurfaceSpecies
from dtcs.common import util
from dtcs.twin import twin_abc
from dtcs.twin.xps import XPSObservable

__author__ = "Andrew Bogdan"
__email__ = "andrewbogdan@lbl.gov"


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

    def __init__(
        self, name: str, binding_energy: float, splitting: float = 1, **kwargs
    ):
        super().__init__(name=name, **kwargs)
        self.binding_energy = binding_energy
        self.splitting = splitting

    def __str__(self):
        if self.splitting == 1:
            return f"{self.name} @ {self.binding_energy} eV"
        else:
            return (
                f"{self.name} @ {self.binding_energy} eV, "
                f"splitting {self.splitting}"
            )

    def __repr__(self):
        if self.splitting == 1:
            return (
                f"{self.__class__.__name__}(name={repr(self.name)}, "
                f"binding_energy={repr(self.binding_energy)})"
            )
        else:
            return (
                f"{self.__class__.__name__}(name={repr(self.name)}, "
                f"binding_energy={repr(self.binding_energy)}, "
                f"splitting={repr(self.splitting)})"
            )


class XPSSpecies(Species):
    """TODO"""

    def __init__(
        self,
        name: str,
        orbs_or_struct: Union[List, Structure] = None,
        *,
        orbitals: List[XPSOrbital] = None,
        structure: Structure = None,
        relax_vis=None,
        xps_vis=None,
        **kwargs,
    ):
        """TODO: Write docstring, check types, make it better, etc."""
        if orbitals:
            if isinstance(orbitals, XPSOrbital):
                orbitals = [orbitals]
            self.orbitals = orbitals
        elif isinstance(orbs_or_struct, (list, tuple)):
            self.orbitals = orbs_or_struct
        elif isinstance(orbs_or_struct, XPSOrbital):
            self.orbitals = [orbs_or_struct]
        elif isinstance(orbs_or_struct, (int, float)):
            self.orbitals = [XPSOrbital(name, binding_energy=orbs_or_struct)]

        if structure:
            self.structure = structure
        elif util.feature_loaded("matproj") and isinstance(orbs_or_struct, Structure):
            self.structure = orbs_or_struct

        if not (hasattr(self, "orbitals") or hasattr(self, "structure")):
            raise TypeError("Either orbitals or structure required.")

        if relax_vis:
            self.relax_vis = relax_vis
        if xps_vis:
            self.xps_vis = xps_vis
        Species.__init__(self, name=name, **kwargs)  # TODO(Andrew): Bodge

    def __str__(self):
        # TODO: will break if there's a structure, not orbitals
        orbitals = [str(orbital) for orbital in self.orbitals]
        if self.color:
            return f"{self.name}, orbitals: {orbitals}, color: {str(self.color)}"
        return f"{self.name}, orbitals: {orbitals}"

    def __repr__(self):
        # TODO: will break if there's a structure, not orbitals
        return (
            f"{self.__class__.__name__}(name={repr(self.name)}, "
            f"orbitals={repr(self.orbitals)}" + f"color={repr(self.color)}"
        )


class XPSSurfaceSpecies(XPSSpecies, SurfaceSpecies):
    def __init__(
        self,
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
        **kwargs,
    ):

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
            **kwargs,
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


class XPSSpec(Spec):
    """A container for a simulated observable of an XPS experiment.

    Attributes:
        autoscale: Defaults to true, decides if it will automatically scale the gaussians and envelope to match the experimental data.
        sim_concs: The concentrations of each species, for the creation of simulated data.
        species_manager: The SpeciesManager in use.
        ignore: Species to ignore during processing.
    """

    def __init__(
        self,
        # XPSObservable args:
        species_manager: SpeciesManager,
        # XPSExperiment args:
        species: Optional[List[sym.Symbol]] = None,
        # The 'cache' for self.resample() to default to.
        x_range: Optional[np.ndarray] = None,
        sim_concs: Optional[Dict[sym.Symbol, float]] = None,
        scale_factor: Optional[float] = None,
        experimental: Optional[pd.Series] = None,
        gas_interval: Optional[Tuple[float, float]] = None,
        contam_spectra: Optional[Dict[sym.Symbol, pd.Series]] = None,
        deconv_species: Optional[List[sym.Symbol]] = None,
        # Flags for how to behave
        autoscale: bool = True,
    ):  # TODO: Add init resample options
        """Initialze the XPSExperiment.

        Pass arguments instead of setting to prevent excessive re-sampling.

        Args: TODO
        """
        super().__init__()

        # --- Parse Arguments ------------------------------------------------
        # Species Manager
        self.species_manager: SpeciesManager = species_manager

        # Require either sim_concs or experimental
        if not (sim_concs or experimental is not None):
            raise ValueError(
                f"{self.__class__.__name__} needs at least"
                f"sim_concs or experimental defined."
            )

        # Make everything match its intended type.
        if x_range is not None:
            x_range = x_range.copy()
        if sim_concs is None:
            sim_concs = {}
        if contam_spectra is None:
            contam_spectra = {}
        if deconv_species is None:
            deconv_species = []
        if experimental is not None:
            experimental = experimental.copy()
        if species is None:
            species = []
            species.extend(sim_concs.keys())
            species.extend(contam_spectra.keys())
            species.extend(deconv_species)
            if not species:
                _logger.echo(
                    f"No species specified to "
                    f"{self.__class__.__name__}(), so all species in "
                    f"the SpeciesManager will be included."
                )
                species = species_manager.species

        # Default scale factor to 1
        # If they pick anything on their own, it disables autoscale.
        if scale_factor is None:
            scale_factor = 1.0
        else:
            autoscale = False

        # Flags
        self.autoscale: bool = autoscale
        # self.species, it will never involve a species outside this list.
        self._species: List[sym.Symbol] = species
        # Internal data, these act like defaults for resampling.
        # If you resample with overwrite=True, they will be overwritten
        # Think of it like a cache.
        self._scale_factor: float = scale_factor
        self._sim_concs: Dict[sym.Symbol, float] = sim_concs
        self._exp_input: Union[pd.Series, None] = experimental
        self._gas_interval: Union[Tuple[float, float], None] = gas_interval
        self._contam_concs: Dict[sym.Symbol, float] = None
        self._contam_spectra: Dict[sym.Symbol, pd.Series] = contam_spectra
        self._deconv_species: List[sym.Symbol] = deconv_species
        self._deconv_concs: Union[Dict[sym.Symbol, float], None] = None

    def simulate(
        self,
        species: Optional[Union[List[sym.Symbol], sym.Symbol]] = None,
        ignore: Optional[Union[List[sym.Symbol], sym.Symbol]] = None,
        x_range: Optional[np.ndarray] = None,
        # Arguments for the simulated data
        simulate: Optional[bool] = None,
        autoscale: Optional[bool] = None,
        scale_factor: Optional[float] = None,
        sim_concs: Optional[dict] = None,
        # Arguments for the experimental data
        experimental: Optional[bool] = None,
        exp_data: Optional[pd.Series] = None,
        gas_phase: Optional[bool] = False,
        gas_interval: Optional[Tuple[float, float]] = None,
        decontaminate: Optional[bool] = None,
        contaminate: Optional[bool] = None,
        contam_spectra: Optional[Dict[sym.Symbol, pd.Series]] = None,
        deconvolute: Optional[bool] = None,
        deconv_species: Optional[List[sym.Symbol]] = None,
        augment=False,
    ) -> XPSObservable:
        """TODO"""
        # --- Parse Arguments ------------------------------------------------
        # Handle: species, ignore
        species = twin_abc._get_species_not_ignored(species, ignore, self.species)

        # --- Simulation-Related ---------------------------------------------
        # Handle: species (Create: orbitals)
        orbitals = {}
        for specie in self.species_manager:
            orbitals[specie.name] = []
            for orbital in specie.orbitals:
                orbitals[specie.name].append(
                    (orbital.binding_energy, orbital.splitting)
                )

        # Handle: sim_concs
        if sim_concs is None:
            sim_concs = self._sim_concs
        # Only take into account species given by _get_species_not_ignored
        sim_species = [specie for specie in sim_concs.keys() if specie in species]

        # Handle: simulate
        # Require that sim_species is nonempty
        if simulate and not sim_species:
            raise ValueError(
                "Cannot simulate without simulated concentrations "
                "defined. Input sim_concs."
            )
        elif simulate is None:
            # Default to whether or not sim_species is empty
            simulate = bool(sim_species)

        # Handle: autoscale, scale_factor
        # scale_factor and autoscale are only used in simulated data.
        if scale_factor is not None and not simulate:
            warnings.warn(
                UserWarning(
                    "Scaling is only used in simulated data;"
                    " user specified scale_factor, but there"
                    " is no simulated data, so it will not"
                    "be used."
                )
            )
        if autoscale and not simulate:
            warnings.warn(
                UserWarning(
                    "Scaling is only used in simulated data;"
                    " user specified autoscale=True, but "
                    "there is no simulated data, so it will "
                    "not be used."
                )
            )
        # Default scale_factor to self._scale_factor/autoscale.
        if scale_factor is None:
            # If autoscale was not given, default to self.autoscale
            if autoscale is None:
                autoscale = self.autoscale
            scale_factor = self._scale_factor  # Won't get used if autoscale.
        elif autoscale:
            raise ValueError("Cannot specify both scale_factor and " "autoscale=True.")

        # --- Experimental-Related -------------------------------------------
        # Handle: exp_data
        if exp_data is None:
            exp_data = self._exp_input  # Might still be None

        # Handle: experimental
        # Experimental requires exp_data
        if experimental and exp_data is None:
            raise ValueError(
                "Experimental data required for experimental-" "related functions."
            )
        elif experimental is None:
            # Default to whether or not there's experimental data.
            experimental = exp_data is not None

        # Handle: gas_interval
        if gas_interval and not experimental:
            warnings.warn(
                UserWarning("Experimental required for gas phase " "evaluation.")
            )
        elif gas_interval is None:
            gas_interval = self._gas_interval

        # Handle: gas_phase
        # gas_phase requires a valid gas_interval
        def valid_gas_interval(interval):
            return len(interval) == 2 and interval[0] <= interval[1]

        if gas_phase and experimental:
            raise ValueError("Cannot evaluate gas phase without experimental.")
        elif gas_phase and not valid_gas_interval(gas_interval):
            raise ValueError(f"Invalid gas interval {gas_interval}")

        # Handle: contam_spectra
        if contam_spectra and not experimental:
            warnings.warn(
                UserWarning("Experimental data required for " "contaminant evaluation.")
            )
        if contam_spectra is None:
            contam_spectra = self._contam_spectra
        # Only take into account species given by _get_species_not_ignored
        contam_species = [
            specie for specie in contam_spectra.keys() if specie in species
        ]  # Everything subsets species

        # Handle: decontaminate, contaminate
        # (de)contaminate requires experimental
        if decontaminate and not experimental:
            raise ValueError("Decontamination requires experimental data.")
        elif contaminate and not simulate:
            raise ValueError("Contamination requires simulated data.")
        # (de)contaminate requires contam_species to be nonempty
        elif decontaminate and not contam_species:
            raise ValueError("Decontamination requries contaminants.")
        elif contaminate and not contam_species:
            raise ValueError("Contamination requires contaminants.")
        # Doing both makes garbage data
        elif decontaminate and contaminate:
            warnings.warn(
                UserWarning(
                    "Both contaminate and decontaminate are "
                    "True, this will yield nonsensical "
                    "data."
                )
            )
        elif contam_species:
            # decontaminate defaults to experimental, contaminate to simulated,
            # and they should not both be True. Decontaminate gets first pick
            # because it's a more common operation.
            if decontaminate is None:
                decontaminate = experimental and not contaminate
            if contaminate is None:
                contaminate = simulate and not decontaminate
        else:
            contaminate = False
            decontaminate = False

        # Handle: deconv_species
        if deconv_species is None:
            deconv_species = self._deconv_species  # TODO: is this dumb?
            # If the default is unspecified, default to species
            if not deconv_species:
                deconv_species = species
        # Only take into account species given by _get_species_not_ignored
        deconv_species = [specie for specie in deconv_species if specie in species]
        # Don't deconvolute species which are being decontaminated out
        deconv_species = [
            specie for specie in deconv_species if specie not in contam_species
        ]
        # Handle: deconvolute
        # deconvolute requires experimental
        if deconvolute and not experimental:
            raise ValueError("Deconvolution requires experimental.")
        # deconvolute requires deconv_species to be nonempty
        elif deconvolute and not deconv_species:
            raise ValueError("Deconvolution requires more than zero species.")
        elif deconvolute is None:
            # Default deconvolute to experimental
            deconvolute = experimental
            if deconvolute:
                deconvolute = _logger.prompt_yn(
                    "Deconvolute experimental?", default=True
                )

        # --- X-Range --------------------------------------------------------
        # Handle: x_range
        # If the user gave an x_range, use it unless there's experimental data
        # and the sizes don't match.
        if x_range is not None and exp_data is not None:
            if x_range.size != exp_data.size:
                raise ValueError(
                    f"Experimental data and supplied x-range have"
                    f"mismatched dimensions {x_range.size}, "
                    f"{exp_data.size}"
                )
        elif exp_data is not None:
            # If the user gave no x_range, default to exp_data.index
            x_range = exp_data.index
        # x_range might still be None, in which case _resample() will scale it
        # automatically.

        # --- Resample and Overwrite -----------------------------------------
        df, scale_factor, deconv_concs, contam_concs = simulate_xps(
            x_range=x_range,
            species=species,
            orbitals=orbitals,
            simulate=simulate,
            sim_concs=sim_concs,
            sim_species=sim_species,
            autoscale=autoscale,
            scale_factor=scale_factor,
            experimental=experimental,
            exp_data=exp_data,
            gas_phase=gas_phase,
            gas_interval=gas_interval,
            decontaminate=decontaminate,
            contaminate=contaminate,
            contam_spectra=contam_spectra,
            contam_species=contam_species,
            deconvolute=deconvolute,
            deconv_species=deconv_species,
            augment=augment,
        )

        return XPSObservable(species_manager=self.species_manager, df=df)

    # TODO(Andrew)
    @classmethod
    def from_cts(cls, cts):
        pass

    # --- Accessors ----------------------------------------------------------

    @property
    def scale_factor(self) -> float:
        """The factor by which the simulated data is currently scaled."""
        return self._scale_factor

    @scale_factor.setter
    def scale_factor(self, scale_factor):
        """The factor by which the simulated data is currently scaled."""
        self._scale_factor = scale_factor

    @property
    def peaks(self) -> Optional[List[float]]:
        peaks = []
        for specie in self.species:
            # TODO(Andrew): This doesn't handle splitting
            bind_eng = self.species_manager[specie].orbitals[0].binding_energy
            peak_index = np.abs(bind_eng - self.x_range).argmin()
            peaks.append(self.x_range[peak_index])
        return peaks

    @property
    def sim_concs(self) -> Optional[Dict[sym.Symbol, float]]:
        """The simulated concentrations of each species."""
        return self._sim_concs

    @property
    def exp_input(self) -> Optional[pd.Series]:
        """The raw experimental data. Might be None.

        Overrides XPSObservable.experimental_raw, as otherwise resampling
        could make it none.
        """
        return self._exp_input

    @property
    def gas_interval(self) -> Optional[Tuple[float, float]]:
        """The interval in which the gas phase peak should be. May be None."""
        return self._gas_interval

    @property
    def contam_concs(self) -> Optional[Dict[sym.Symbol, float]]:
        """TODO"""
        return self._contam_concs

    @property
    def contam_spectra(self) -> Optional[Dict[sym.Symbol, pd.Series]]:
        """The spectra of the contaminants used to decontaminate the data.

        Won't usually match x_range.
        """
        return self._contam_spectra

    @property
    def deconv_species(self) -> Optional[List[sym.Symbol]]:
        """TODO"""
        return self._deconv_species

    @property
    def deconv_concs(self) -> Optional[Dict[sym.Symbol, float]]:
        """The concentrations of each species guessed by the deconvolution."""
        return self._deconv_concs

    @property
    def species(self) -> List[sym.Symbol]:
        """All the species in the experiment."""
        return self._species

    # --- Utility -------------------------------------------------------------

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict()
        d["sim_concs"] = {str(symbol): conc for symbol, conc in self._sim_concs.items()}
        d["species_manager"] = self.species_manager.as_dict()
        d["autoscale"] = self.autoscale
        d["experimental"] = (
            self._exp_input.to_json() if self._exp_input is not None else None
        )
        d["gas_interval"] = self._gas_interval
        d["scale_factor"] = self._scale_factor
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        decode = monty.json.MontyDecoder().process_decoded
        d["sim_concs"] = {
            sym.Symbol(name): conc for name, conc in d["sim_concs"].items()
        }
        d["species_manager"] = decode(d["species_manager"])
        if d["experimental"] is not None:
            d["experimental"] = pd.read_json(
                d["experimental"], typ="series", convert_axes=False
            )
            d["experimental"].index = d["experimental"].index.map(float)
        return cls(**d)


# TODO: For backwards compatibility
Orbital = XPSOrbital
