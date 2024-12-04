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
from dtcs.sim.ir import simulate_ir

from dtcs.spec.spec_abc import Spec
from dtcs.spec.species import Species, SpeciesManager
from dtcs.spec.crn.surface.species import SurfaceSpecies
from dtcs.common import util
from dtcs.twin import twin_abc
from dtcs.twin.xps import XPSObservable

from dtcs.spec.ir import IRSpecies, IRSignature


class IRSurfaceSpecies(IRSpecies, SurfaceSpecies):
    def __init__(
        self,
        name: str,
        ir_signature: Tuple[List[IRSignature] | IRSignature],
        site=None,
        size=None,
        is_gas=False,
        is_visible=True,
        **kwargs,
    ):

        self.surface_species = SurfaceSpecies.__init__(
            self,
            name=name,
            site=site,
            size=size,
            is_gas=is_gas,
            is_visible=is_visible,
        )

        self.ir_species = IRSpecies.__init__(
            self, name=name, ir_signature=ir_signature, is_visible=is_visible
        )


class IRSpec(Spec):
    """A container for a simulated observable of an IR experiment.

    Attributes:
        autoscale: Defaults to true, decides if it will automatically scale the gaussians and envelope to match the experimental data.
        sim_concs: The concentrations of each species, for the creation of simulated data.
        species_manager: The SpeciesManager in use.
        ignore: Species to ignore during processing.
    """

    def __init__(
        self,
        species_manager: SpeciesManager,
        species: Optional[List[sym.Symbol]] = None,
        x_range: Optional[np.ndarray] = None,
        sim_concs: Optional[Dict[sym.Symbol, float]] = None,
        scale_factor: Optional[float] = None,
        experimental: Optional[pd.Series] = None,
        gas_interval: Optional[Tuple[float, float]] = None,
        deconv_species: Optional[List[sym.Symbol]] = None,
        autoscale: bool = True,
    ):
        """Initialze the IRExperiment.
        Pass arguments instead of setting to prevent excessive re-sampling.
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
        if deconv_species is None:
            deconv_species = []
        if experimental is not None:
            experimental = experimental.copy()
        if species is None:
            species = []
            species.extend(sim_concs.keys())
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

        # Handle: deconv_species
        if deconv_species is None:
            deconv_species = self._deconv_species  # TODO: is this dumb?
            # If the default is unspecified, default to species
            if not deconv_species:
                deconv_species = species
        # Only take into account species given by _get_species_not_ignored
        deconv_species = [specie for specie in deconv_species if specie in species]

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
        df, scale_factor, deconv_concs = simulate_ir(
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
