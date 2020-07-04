"""Utilities for manipulating XPS data.

Exports:
    XPSExperiment: a container for simulated and experimental xps observables.

Example:
    Say you have a normal setup with a species manager, species x, y, and some
    experimental data in a pd.Series. The following is preferred:

    xps = XPSExperiment({x: 2.1, y: 3.4}, sm, experimental=experimental_data,
                        gas_interval=(-1, 0), scale_factor=0.1)
    xps.plot()

    but is equivalent to the following, which resamples with each assignment:

    xps = XPSExperiment({x: 2.1, y: 3.4}, sm)
    xps.experimental = experimental_data
    xps.gas_interval = (-1, 0)
    xps.scale_factor = 0.1
    xps.plot()
"""

import copy
from typing import Dict, List, Optional, Tuple, Union
import warnings

from matplotlib import pyplot as plt
import monty.json
import numpy as np
import pandas as pd
from scipy import integrate
from scipy import optimize
from scipy import stats
from sklearn import metrics
import sympy as sym

from lblcrn import bulk_crn
from lblcrn.experiments import experiment
from lblcrn.experiments import time_series
from lblcrn.crn_sym import reaction
from lblcrn.crn_sym import species
from lblcrn.common import util
from lblcrn import _echo


class XPSObservable:
    """TODO"""

    # Settings for the column names
    _SIMULATED = 'simulated'
    _EXPERIMENTAL = 'experimental'
    _CONTAMINANTS = 'contaminants'
    _DECONVOLUTED = 'deconvoluted'
    _SIM_ENV = (_SIMULATED, 'envelope')
    _CLEAN_EXP = (_EXPERIMENTAL, 'clean')
    _RAW_EXP = (_EXPERIMENTAL, 'raw')
    _GAS_PHASE = (_EXPERIMENTAL, 'gas_phase')
    _DECONV_ENV = (_DECONVOLUTED, 'envelope')

    def __init__(self, df: pd.DataFrame,
                 species_manager: species.SpeciesManager,
                 title: str = ''):
        self.df = df
        self.species_manager = species_manager
        self.title = title

    @property
    def x_range(self) -> np.ndarray:
        """The x-values, energies, on which there is data."""
        return np.asarray(self.df.index)

    # TODO: Perhaps add an 'envelopes' property. It would be useful for
    #  error analysis and changing the opacity/line-width of all the
    #  envelopes at once.

    @property
    def gaussians(self) -> Optional[pd.DataFrame]:
        """The gaussians of the XPS observable."""
        # Columns are a pair (simulated/experimental, str/species)
        # If they're a species, we want them.
        gauss_cols = []
        if self.gas_phase is not None:
            gauss_cols.append(self._GAS_PHASE)
        for col in self.df.columns:
            if isinstance(col[1], sym.Symbol):
                gauss_cols.append(col)
        if gauss_cols:
            return self.df[gauss_cols]
        else:
            return None

    @property
    def simulated(self) -> Optional[pd.DataFrame]:
        """The simulated data. Might be None."""
        if self._SIMULATED in self.df:
            return self.df[self._SIMULATED]
        return None

    @property
    def sim_gaussians(self) -> Optional[pd.DataFrame]:
        """The simulated gaussians of the XPS observable."""
        # Gaussians have a species as their column name
        gauss_cols = []
        for col in self.simulated:
            if isinstance(col, sym.Symbol):
                gauss_cols.append(col)
        if gauss_cols:
            return self.simulated[gauss_cols]
        else:
            return None

    @property
    def sim_envelope(self) -> Optional[pd.Series]:
        """The simulated envelope, the sum of all species' gaussians."""
        if self._SIM_ENV in self.df:
            return self.df[self._SIM_ENV]
        return None

    @property
    def experimental(self) -> Optional[pd.DataFrame]:
        """The experimental data. Might be None."""
        if self._EXPERIMENTAL in self.df:
            return self.df[self._EXPERIMENTAL]
        return None

    @property
    def exp_clean(self) -> Optional[pd.Series]:
        """The experimental data without the gas phase or contaminants."""
        if self._CLEAN_EXP in self.df:
            return self.df[self._CLEAN_EXP]
        return None

    @property
    def exp_raw(self) -> Optional[pd.Series]:
        """The raw experimental data. Might be None."""
        if self._RAW_EXP in self.df:
            return self.df[self._RAW_EXP]
        return None

    @property
    def gas_phase(self) -> Optional[pd.Series]:
        """The gas phase part of the spectrum. Might be None."""
        if self._GAS_PHASE in self.df:
            return self.df[self._GAS_PHASE]
        return None

    @property
    def contaminants(self) -> Optional[pd.DataFrame]:
        """The contaminants' contribution to the spectrum."""
        if self._CONTAMINANTS in self.df:
            return self.df[self._CONTAMINANTS]
        return None

    @property
    def deconvoluted(self) -> Optional[pd.DataFrame]:
        """All of the deconvolution data."""
        if self._DECONVOLUTED in self.df:
            return self.df[self._DECONVOLUTED]
        return None

    @property
    def deconv_gaussians(self) -> Optional[pd.DataFrame]:
        """The deconvoluted gaussians of the XPS observable."""
        # Gaussians have a species as their column name
        gauss_cols = []
        for col in self.deconvoluted:
            if isinstance(col, sym.Symbol):
                gauss_cols.append(col)
        if gauss_cols:
            return self.deconvoluted[gauss_cols]
        else:
            return None

    @property
    def deconv_envelope(self) -> Optional[pd.Series]:
        """The envelope of the deconvoluted experimental data."""
        if self._DECONV_ENV in self.df:
            return self.df[self._DECONV_ENV]
        return None

    def plot(self, ax: plt.Axes = None,
             only: bool = False,
             species: Optional[Union[List[sym.Symbol], sym.Symbol]] = None,
             ignore: Optional[Union[List[sym.Symbol], sym.Symbol]] = None,
             peak_lines: Optional[Union[dict, bool]] = None,
             gaussians: Optional[Union[dict, bool]] = None,
             simulated: Optional[Union[dict, bool]] = None,
             sim_envelope: Optional[Union[dict, bool]] = None,
             sim_gaussians: Optional[Union[dict, bool]] = None,
             experimental: Optional[Union[dict, bool]] = None,
             exp_clean: Optional[Union[dict, bool]] = None,
             exp_raw: Optional[Union[dict, bool]] = None,
             gas_phase: Optional[Union[dict, bool]] = None,
             contaminants: Optional[Union[dict, bool]] = None,
             deconvoluted: Optional[Union[dict, bool]] = None,
             deconv_envelope: Optional[Union[dict, bool]] = None,
             deconv_gaussians: Optional[Union[dict, bool]] = None,
             ) -> plt.Axes:
        """TODO"""
        # --- Parse Arguments ------------------------------------------------
        if ax is None:
            ax = plt.gca()

        species = experiment._get_species_not_ignored(
            species,
            ignore,
            self.species_manager.species,
        )

        # Some switches count for more than one thing.
        #  We ignore these switches if they aren't explicitly specified.
        # I don't care if this goes against PEP 8 this would be a mess!
        if gaussians is not None:
            default = gaussians or isinstance(gaussians, dict)
            if sim_gaussians is None: sim_gaussians = default
            if gas_phase is None: gas_phase = default
            if contaminants is None: contaminants = default
            if deconv_gaussians is None: deconv_gaussians = default
        if simulated is not None:
            default = simulated or isinstance(simulated, dict)
            if sim_gaussians is None: sim_gaussians = default
            if sim_envelope is None: sim_envelope = default
        if experimental is not None:
            default = experimental or isinstance(experimental, dict)
            if exp_raw is None: exp_raw = default
            if exp_clean is None: exp_clean = default
            if gas_phase is None: gas_phase = bool(experimental)
        if deconvoluted is not None:
            default = deconvoluted or isinstance(deconvoluted, dict)
            if deconv_gaussians is None: deconv_gaussians = default
            if deconv_envelope is None: deconv_envelope = default

        if peak_lines is None: peak_lines = not only
        if sim_envelope is None: sim_envelope = not only
        if sim_gaussians is None: sim_gaussians = not only
        if exp_clean is None: exp_clean = not only
        if exp_raw is None: exp_raw = not only
        if gas_phase is None: gas_phase = not only
        if contaminants is None: contaminants = not only
        if deconv_envelope is None: deconv_envelope = not only
        if deconv_gaussians is None: deconv_gaussians = not only

        # Get the gaussians we're going to be plotting
        sim_gauss_species = []
        deconv_gauss_species = []
        contam_gauss_species = []
        species_to_plot = []
        if sim_gaussians and self.sim_gaussians is not None:
            sim_gauss_species.extend(col for col in self.sim_gaussians
                                  if col in species)
        if deconv_gaussians and self.deconv_gaussians is not None:
            deconv_gauss_species.extend(col for col in self.deconv_gaussians
                                     if col in species)
        if contaminants and self.contaminants is not None:
            contam_gauss_species.extend(specie for specie in self.contaminants
                                     if specie in species)
        species_to_plot.extend(sim_gauss_species)
        species_to_plot.extend(deconv_gauss_species)
        species_to_plot.extend(contam_gauss_species)

        # --- Plotting Arguments ---------------------------------------------
        # When plotting simulated/deconvoluted species, we're going to have
        #  gaussians on top of other gaussians which might be the same color.
        #  Hence, we have different arguments for gaussians which are plotted
        #  first, second, etc.
        ordered_gauss_args = [
            dict(fill=True, ),
            dict(fill=False, brightness=0.8, ),
            dict(fill=False, brightness=0.6, ),
            dict(fill=False, brightness=0.4, ),
        ]

        # Arguments for groups
        gaussians_args = dict()
        simulated_args = dict()
        experimental_args = dict()
        deconvoluted_args = dict()

        # Args: peak_lines
        peak_line_args = dict(
            ymin=0.05,
            ymax=1,
            linestyle='dashed',
            linewidth=1
        )
        # Args: sim_envelope
        envelope_args = copy.copy(simulated_args)
        envelope_args.update(
            label='Sim. Envelope',
            color='#663d3d',
            linestyle='--',
            linewidth=3,
        )
        # Args: exp_clean
        exp_clean_args = copy.copy(experimental_args)
        exp_clean_args.update(
            label='Clean Experimental',
            color='0',
            linewidth=3,
        )
        # Args: exp_raw
        exp_raw_args = copy.copy(experimental_args)
        exp_raw_args.update(
            label='Raw Experimental',
            color='0.3',
            linestyle='-.',
            linewidth=3,
        )
        # Args: gas_phase
        gas_phase_args = copy.copy(gaussians_args)
        gas_phase_args.update(experimental_args)
        gas_phase_args.update(
            label='Gas Phase',
            color='#666666',
            fill=False,
            hatch='++'
        )
        # Args: deconv_envelope
        deconv_env_args = copy.copy(deconvoluted_args)
        deconv_env_args.update(
            label='Deconv. Envelope',
            color='#3d663d',
            linestyle=':',
            linewidth=3,
        )

        # These should be in plotting order, or else ordered_gauss_args
        #  won't work.
        # Args: sim_gaussians
        sim_gauss_args = copy.copy(gaussians_args)
        if sim_gauss_species:
            sim_gauss_args.update(ordered_gauss_args.pop(0))
        sim_gauss_args.update(simulated_args)
        sim_gauss_args.update(
            hatch='xxx',
        )
        # Args: deconv_gaussians
        deconv_gauss_args = copy.copy(gaussians_args)
        if deconv_gauss_species:
            deconv_gauss_args.update(ordered_gauss_args.pop(0))
        deconv_gauss_args.update(deconvoluted_args)
        deconv_gauss_args.update(
            hatch='///',
        )
        # Args: contaminants
        contam_args = copy.copy(gaussians_args)
        if contam_gauss_species:
            contam_args.update(ordered_gauss_args.pop(0))
        contam_args.update(
            hatch='...',
        )

        # If the user gave a dictionary, use those.
        if isinstance(gaussians, dict):
            gas_phase_args.update(gaussians)
            sim_gauss_args.update(gaussians)
            deconv_gauss_args.update(gaussians)
            contam_args.update(gaussians)
        if isinstance(simulated, dict):
            sim_gauss_args.update(simulated)
            envelope_args.update(simulated)
        if isinstance(experimental, dict):
            gas_phase_args.update(experimental)
            exp_raw_args.update(experimental)
            exp_clean_args.update(experimental)
        if isinstance(deconvoluted, dict):
            deconv_gauss_args.update(deconvoluted)
            deconv_env_args.update(deconvoluted)

        for inputted, args_dict in [(peak_lines, peak_line_args),
                                    (sim_gaussians, sim_gauss_args),
                                    (sim_envelope, envelope_args),
                                    (exp_raw, exp_raw_args),
                                    (exp_clean, exp_clean_args),
                                    (gas_phase, gas_phase_args),
                                    (contaminants, contam_args),
                                    (deconv_gaussians, deconv_gauss_args),
                                    (deconv_envelope, deconv_env_args)]:
            if isinstance(inputted, dict):
                args_dict.update(inputted)

        # --- Plot -----------------------------------------------------------
        # Function to retrieve with arguments not to be passed to pyplot,
        #  you must let it modify the dict you give it, or else it'll pass
        #  bad arguments to pyplot.
        def pop_special_args(args_dict):
            brightness = 1.0
            if 'brightness' in args_dict:
                brightness = args_dict.pop('brightness')
            return brightness

        if peak_lines and species_to_plot:
            brightness = pop_special_args(peak_line_args)
            for specie in species_to_plot:
                for orbital in self.species_manager[specie].orbitals:
                    if orbital.binding_energy not in self.x_range:
                        continue
                    color = self.species_manager.color(specie)
                    color = util.scale_color_brightness(color, brightness)

                    ax.axvline(x=orbital.binding_energy,
                               color=color,
                               **peak_line_args)

        if gas_phase and self.gas_phase is not None:
            _ = pop_special_args(gas_phase_args)
            ax.fill(self.x_range, self.gas_phase, **gas_phase_args)

        if sim_gaussians and sim_gauss_species:
            brightness = pop_special_args(sim_gauss_args)

            # Sort so that the shorter peaks are first
            def get_max_of_sim_column(specie):
                return max(self.sim_gaussians[specie])
            
            for specie in sorted(sim_gauss_species,
                                 key=get_max_of_sim_column,
                                 reverse=True):
                color = self.species_manager.color(specie)
                color = util.scale_color_brightness(color, brightness)

                ax.fill(self.x_range, self.sim_gaussians[specie],
                        label=f'Sim. {specie}',
                        color=color,
                        **sim_gauss_args)

        if deconv_gaussians and deconv_gauss_species:
            brightness = pop_special_args(deconv_gauss_args)

            # Sort so that the shorter peaks are first
            def get_max_of_deconv_column(specie):
                return max(self.deconv_gaussians[specie])
            
            for specie in sorted(deconv_gauss_species,
                                 key=get_max_of_deconv_column,
                                 reverse=True):
                color = self.species_manager.color(specie)
                color = util.scale_color_brightness(color, brightness)

                ax.fill(self.x_range, self.deconv_gaussians[specie],
                        label=f'Deconv. {specie}',
                        color=color,
                        **deconv_gauss_args)

        if contaminants and contam_gauss_species:
            brightness = pop_special_args(contam_args)

            # Sort so that the shorter peaks are first
            def get_max_of_contam_column(specie):
                return max(self.contaminants[specie])

            for specie in sorted(contam_gauss_species,
                                 key=get_max_of_contam_column,
                                 reverse=True):
                color = self.species_manager.color(specie)
                color = util.scale_color_brightness(color, brightness)

                ax.fill(self.x_range, self.contaminants[specie],
                        label=f'Contam. {specie}',
                        color=color,
                        **contam_args)

        if deconv_envelope and self.deconv_envelope is not None:
            _ = pop_special_args(deconv_env_args)
            ax.plot(self.x_range, self.deconv_envelope, **deconv_env_args)

        if sim_envelope and self.sim_envelope is not None:
            _ = pop_special_args(envelope_args)
            ax.plot(self.x_range, self.sim_envelope, **envelope_args)

        if exp_raw and self.exp_raw is not None:
            _ = pop_special_args(exp_raw_args)
            ax.plot(self.x_range, self.exp_raw, **exp_raw_args)

        if exp_clean and self.exp_clean is not None:
            _ = pop_special_args(exp_clean_args)
            ax.plot(self.x_range, self.exp_clean, **exp_clean_args)

        ax.legend()
        # TODO: This overrides whatever title they set for the axis
        #  beforehand. Make it an option which XPSExperiment/CRNTimeSeries
        #  passes itself?
        ax.set_title(self.title)
        ax.invert_xaxis()  # XPS Plots are backwards
        return ax

    # --- Data Analysis ------------------------------------------------------

    def rmse(self):
        return np.sqrt(metrics.mean_squared_error(self.exp_raw,
                                                  self.sim_envelope))

    def mae(self):
        return metrics.mean_absolute_error(self.exp_raw,
                                           self.sim_envelope)

    def integral_diff_outside_experimental(self):
        return integrate.trapz(self.sim_envelope - self.exp_raw,
                               self.x_range)

    def integral_diff_inside_experimental(self):
        return integrate.trapz(self.exp_raw - self.sim_envelope,
                               self.x_range)

    def integral_diff_between(self):
        return integrate.trapz(np.absolute(self.exp_raw -
                                           self.sim_envelope), self.x_range)

    def envelope_integral(self):
        return integrate.trapz(self.sim_envelope, self.x_range)


class XPSExperiment(experiment.Experiment, XPSObservable):
    """A container for a simulated observable of an XPS experiment.

    Attributes:
        autoresample: Defaults to true, decides if it resamples on edits.
        autoscale: Defaults to true, decides if it will automatically scale
            the gaussians and envelope to match the experimental data.
        sim_concs: The concentrations of each species, for the creation of
            simulated data.
        species_manager: The SpeciesManager in use.
    """

    _PLOT_MARGIN = 5
    _PLOT_RESOLUTION = 0.001
    _SIGMA = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)

    def __init__(self,
                 # XPSObservable args:
                 species_manager: species.SpeciesManager,
                 title: Optional[str] = '',
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
                 autoresample: bool = True,
                 autoscale: bool = True, ):  # TODO: Add init resample options
        """Initialze the XPSExperiment.

        Pass arguments instead of setting to prevent excessive re-sampling.

        Args: TODO
        """
        super().__init__()
        XPSObservable.__init__(self, species_manager=species_manager,
                               df=None, title=title)

        # --- Parse Arguments ------------------------------------------------
        # Require either sim_concs or experimental
        if not (sim_concs or experimental is not None):
            raise ValueError(f'{self.__class__.__name__} needs at least'
                             f'species_concs or experimental defined.')

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
                _echo.echo(f'No species specified to '
                           f'{self.__class__.__name__}(), so all species in '
                           f'the SpeciesManager will be included.')
                species = species_manager.species

        # Default scale factor to 1
        # If they pick anything on their own, it disables autoscale.
        if scale_factor is None:
            scale_factor = 1.0
        else:
            autoscale = False

        # Flags
        self.autoresample: bool = autoresample
        self.autoscale: bool = autoscale
        # self.species, it will never involve a species outside this list.
        self._species: List[sym.Symbol] = species
        # Internal data, these act like defaults for resampling.
        # If you resample with overwrite=True, they will be overwritten
        # Think of it like a cache.
        self._x_range: Union[np.ndarray, None] = x_range
        self._sim_concs: Dict[sym.Symbol, float] = sim_concs
        self._scale_factor: float = scale_factor
        self._exp_data: Union[pd.Series, None] = experimental
        self._gas_interval: Union[Tuple[float, float], None] = gas_interval
        self._contam_spectra: Dict[sym.Symbol, pd.Series] = contam_spectra
        self._deconv_species: List[sym.Symbol] = deconv_species
        self._deconv_concs: Union[Dict[sym.Symbol, float], None] = None

        # Resample if autoresample is True.
        self._autoresample()

    # --- Accessors ----------------------------------------------------------

    @property
    def exp_raw(self) -> Optional[pd.Series]:
        """The raw experimental data. Might be None.

        Overrides XPSObservable.experimental_raw, as otherwise resampling
        could make it none.
        """
        return self._exp_data

    @property
    def gas_interval(self) -> Optional[Tuple[float, float]]:
        """The interval in which the gas phase peak should be. May be None."""
        return self._gas_interval

    @property
    def scale_factor(self) -> float:
        """The factor by which the simulated data is currently scaled."""
        return self._scale_factor

    @property
    def simulated_concs(self) -> Optional[Dict[sym.Symbol, float]]:
        """The simulated concentrations of each species."""
        return self._sim_concs

    @property
    def deconvoluted_concs(self) -> Optional[Dict[sym.Symbol, float]]:
        """The concentrations of each species guessed by the deconvolution."""
        return self._deconv_concs

    @property
    def contaminant_spectra(self) -> Optional[Dict[sym.Symbol, pd.Series]]:
        """The spectra of the contaminants used to decontaminate the data.

        Won't usually match x_range.
        """
        return self._contam_spectra

    @property
    def species(self) -> List[sym.Symbol]:
        """All the species in the experiment."""
        return self._species

    # --- Mutators -----------------------------------------------------------
    # If autoresample is on, these will prompt an overwriting resample.

    def include(self, *species: sym.Symbol):
        """Add the species to the experiment, if they aren't there already.

        It won't simulate data for them (unless you also edit species_concs),
        but it will involve the species in deconvolution and the like.

        Prompts an autoresample.
        """
        for specie in species:
            if specie not in self._species:
                self._species.append(specie)
        self._autoresample()

    def drop(self, *species: sym.Symbol):
        """Remove the species from the experiment.

        It won't deconvolute with or simulate data from them.

        Prompts an autoresample.
        """
        for specie in species:
            self._species.remove(specie)
        self._autoresample()

    def _autoresample(self):
        """Does an overwriting resample, if necessary."""
        if self.autoresample:
            _echo.echo('Auto-resampling data...')
            self.resample(overwrite=True)

    # --- Calculations -------------------------------------------------------
    # Note that functions in this section will only modify state if
    # overwwrite=True.

    def resample(
            self,
            overwrite: bool = True,
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
            gas_phase: Optional[bool] = None,
            gas_interval: Optional[Tuple[float, float]] = None,
            decontaminate: Optional[bool] = None,
            contaminate: Optional[bool] = None,
            contam_spectra: Optional[Dict[sym.Symbol, pd.Series]] = None,
            deconvolute: Optional[bool] = None,
            deconv_species: Optional[List[sym.Symbol]] = None,
    ) -> XPSObservable:
        """Re-calculate everything.

        Args:
            overwrite:
            species:
            ignore:

            --- Data Inputs ---
            x_range:
            scale_factor:
            sim_concs:
            exp_data:
            gas_interval:
            contam_spectra:
            deconv_species:

            --- Flags ---
            simulate:
            autoscale:
            experimental:
            gas_phase:
            decontaminate:
            contaminate:
            deconvolute:

        TODO: Cache all the data inputs

        Returns:

        """
        # --- Parse Arguments ------------------------------------------------
        # Handle: species, ignore
        species = experiment._get_species_not_ignored(species, ignore,
                                                      self.species)

        # --- Simulation-Related ---------------------------------------------
        # Handle: sim_concs
        if sim_concs is None:
            sim_concs = self._sim_concs
        # Only take into account species given by _get_species_not_ignored
        sim_species = [specie for specie in sim_concs.keys()
                       if specie in species]

        # Handle: simulate
        # Require that sim_species is nonempty
        if simulate and not sim_species:
            raise ValueError('Cannot simulate without simulated concentrations '
                             'defined. Input sim_concs.')
        elif simulate is None:
            # Default to whether or not sim_species is empty
            simulate = bool(sim_species)

        # Handle: autoscale, scale_factor
        # scale_factor and autoscale are only used in simulated data.
        if scale_factor is not None and not simulate:
            warnings.warn(UserWarning('Scaling is only used in simulated data;'
                                      ' user specified scale_factor, but there'
                                      ' is no simulated data, so it will not'
                                      'be used.'))
        if autoscale and not simulate:
            warnings.warn(UserWarning('Scaling is only used in simulated data;'
                                      ' user specified autoscale=True, but '
                                      'there is no simulated data, so it will '
                                      'not be used.'))
        # Default scale_factor to self._scale_factor/autoscale.
        if scale_factor is None:
            # If autoscale was not given, default to self.autoscale
            if autoscale is None:
                autoscale = self.autoscale
            scale_factor = self._scale_factor  # Won't get used if autoscale.
        elif autoscale:
            raise ValueError('Cannot specify both scale_factor and '
                             'autoscale=True.')

        # --- Experimental-Related -------------------------------------------
        # Handle: exp_data
        if exp_data is None:
            exp_data = self._exp_data  # Might still be None

        # Handle: experimental
        # Experimental requires exp_data
        if experimental and exp_data is None:
            raise ValueError('Experimental data required for experimental-'
                             'related functions.')
        elif experimental is None:
            # Default to whether or not there's experimental data.
            experimental = (exp_data is not None)

        # Handle: gas_interval
        if gas_interval and not experimental:
            warnings.warn(UserWarning('Experimental required for gas phase '
                                      'evaluation.'))
        elif gas_interval is None:
            gas_interval = self._gas_interval

        # Handle: gas_phase
        # gas_phase requires a valid gas_interval
        def valid_gas_interval(interval):
            return len(interval) == 2 and interval[0] <= interval[1]

        if gas_phase and experimental:
            raise ValueError('Cannot evaluate gas phase without experimental.')
        elif gas_phase and not valid_gas_interval(gas_interval):
            raise ValueError(f'Invalid gas interval {gas_interval}')
        elif gas_phase is None:
            gas_phase = experimental

        # Handle: contam_spectra
        if contam_spectra and not experimental:
            warnings.warn(UserWarning('Experimental data required for '
                                      'contaminant evaluation.'))
        if contam_spectra is None:
            contam_spectra = self._contam_spectra
        # Only take into account species given by _get_species_not_ignored
        contam_species = [specie for specie in contam_spectra.keys()
                          if specie in species]  # Everything subsets species

        # Handle: decontaminate, contaminate
        # (de)contaminate requires experimental
        if decontaminate and not experimental:
            raise ValueError('Decontamination requires experimental data.')
        elif contaminate and not simulate:
            raise ValueError('Contamination requires simulated data.')
        # (de)contaminate requires contam_species to be nonempty
        elif decontaminate and not contam_species:
            raise ValueError('Decontamination requries contaminants.')
        elif contaminate and not contam_species:
            raise ValueError('Contamination requires contaminants.')
        # Doing both makes garbage data
        elif decontaminate and contaminate:
            warnings.warn(UserWarning('Both contaminate and decontaminate are '
                                      'True, this will yield nonsensical '
                                      'data.'))
        else:
            # decontaminate defaults to experimental, contaminate to simulated,
            # and they should not both be True. Decontaminate gets first pick
            # because it's a more common operation.
            if decontaminate is None:
                decontaminate = experimental and not contaminate
            if contaminate is None:
                contaminate = simulate and not decontaminate

        # Handle: deconv_species
        if deconv_species is None:
            deconv_species = self._deconv_species  # TODO: is this dumb?
            # If the default is unspecified, default to species
            if not deconv_species:
                deconv_species = species
        # Only take into account species given by _get_species_not_ignored
        deconv_species = [specie for specie in deconv_species
                          if specie in species]
        # Don't deconvolute species which are being decontaminated out
        deconv_species = [specie for specie in deconv_species
                          if specie not in contam_species]

        # Handle: deconvolute
        # deconvolute requires experimental
        if deconvolute and not experimental:
            raise ValueError('Deconvolution requires experimental.')
        # deconvolute requires deconv_species to be nonempty
        elif deconvolute and not deconv_species:
            raise ValueError('Deconvolution requires more than zero species.')
        elif deconvolute is None:
            # Default deconvolute to experimental
            deconvolute = experimental
            if deconvolute:
                deconvolute = _echo.prompt_yn('Deconvolute experimental?',
                                              default=True)

        # --- X-Range --------------------------------------------------------
        # Handle: x_range
        # If the user gave an x_range, use it unless there's experimental data
        # and the sizes don't match.
        if x_range is not None and exp_data is not None:
            if x_range.size != exp_data.size:
                raise ValueError(f'Experimental data and supplied x-range have'
                                 f'mismatched dimensions {x_range.size}, '
                                 f'{exp_data.size}')
        elif exp_data is not None:
            # If the user gave no x_range, default to exp_data.index
            x_range = exp_data.index
        elif x_range is None:
            # This might be None, in which case _resample() will scale it
            # automatically. Note that _x_range will never be overridden by
            # the automatic x_ranges, only by the user and exp_data.
            x_range = self._x_range

        # --- Resample and Overwrite -----------------------------------------
        df, scale_factor, deconv_concs, contam_concs = self._resample(
            x_range=x_range,
            species=species,
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
        )

        if overwrite:
            self._x_range = x_range
            self._sim_concs = sim_concs
            self._scale_factor = scale_factor
            self._exp_data = exp_data
            self._gas_interval = gas_interval
            self._contam_spectra = contam_spectra
            self._deconv_concs = deconv_concs
            self.df = df
        return XPSObservable(species_manager=self.species_manager,
                             df=df,
                             title=self.title)

    def _resample(
            self,
            # Input Data
            x_range: Union[np.ndarray, None],
            species: List[sym.Symbol],
            sim_concs: Dict[sym.Symbol, float],
            sim_species: List[sym.Symbol],
            scale_factor: float,
            exp_data: pd.Series,
            gas_interval: Tuple[float, float],
            contam_spectra: Dict[sym.Symbol, pd.Series],
            contam_species: List[sym.Symbol],
            deconv_species: List[sym.Symbol],
            # Flags
            simulate: bool,
            autoscale: bool,
            experimental: bool,
            gas_phase: bool,
            decontaminate: bool,
            contaminate: bool,
            deconvolute: bool,
    ) -> Tuple[pd.DataFrame,
               float,
               Union[Dict[sym.Symbol, float], None],
               Union[Dict[sym.Symbol, float], None]]:
        """TODO"""
        # --- Invariants -----------------------------------------------------
        assert not (simulate and not sim_concs)
        assert not (simulate and not sim_species)
        assert not (experimental and exp_data is None)
        assert not (gas_phase and not experimental)
        assert not (gas_phase and (len(gas_interval) > 2 or
                                   gas_interval[0] > gas_interval[1]))
        assert not (decontaminate and not experimental)
        assert not (contaminate and not simulate)
        assert not ((decontaminate or contaminate) and not contam_spectra)
        assert not ((decontaminate or contaminate) and not contam_species)
        assert not (contaminate and decontaminate)
        assert not (deconvolute and not experimental)
        assert not (deconvolute and not deconv_species)
        assert scale_factor != 0

        # --- Variables for Later --------------------------------------------
        envelope: Union[pd.Series, None] = None
        deconv_concs: Union[Dict[sym.Symbol, float], None] = None
        contam_concs: Union[Dict[sym.Symbol, float], None] = None
        sim_cols: List[Tuple[str, sym.Symbol]] = []
        deconv_cols: List[Tuple[str, sym.Symbol]] = []
        contam_cols: List[Tuple[str, sym.Symbol]] = []

        # --- Make X-Range ---------------------------------------------------
        # Find the x_range based on what species we need to see.
        if x_range is None:
            species_to_plot = []
            if simulate:
                species_to_plot.extend(sim_species)
            # Don't bother with deconv_/contam_species, because if they're
            # defined, then the x_range will be exp_data.index anyway.
            x_range = self._get_x_range(species=species_to_plot)
            _echo.echo(f'Using automatically-generated x-range '
                       f'[{x_range[0]}, ..., {x_range[-1]}]...')

        # --- Make Empty DataFrame -------------------------------------------
        # I add all of the columns early for sorting reasons;
        # comparisons between sym.Symbol and str are not defined, so you can't
        # easily sort the columns.
        columns = []
        if simulate:
            columns.append(self._SIM_ENV)
            sim_cols = [(self._SIMULATED, specie) for specie in sim_species]
            columns.extend(sim_cols)
        if experimental:
            columns.append(self._RAW_EXP)
            # Gas phase and decontaminate will clean up the raw experimental.
        if gas_interval or decontaminate:
            columns.append(self._CLEAN_EXP)
        if gas_phase:
            columns.append(self._GAS_PHASE)
        if decontaminate or contaminate:
            contam_cols = [(self._CONTAMINANTS, specie)
                           for specie in contam_species]
        if deconvolute:
            columns.append(self._DECONV_ENV)
            deconv_cols = [(self._DECONVOLUTED, specie)
                           for specie in deconv_species]
            columns.extend(deconv_cols)

        row_index = pd.Index(x_range, name='eV')
        col_index = pd.MultiIndex.from_tuples(columns)
        df = pd.DataFrame(data=0, index=row_index, columns=col_index)

        # --- Calculations ---------------------------------------------------
        if simulate:
            if _echo.do_echo:
                sim_concs_used = {specie: conc for specie, conc in
                                  sim_concs.items() if specie in sim_species}
                _echo.echo(f'Simulating XPS experiment for species with '
                           f'concentrations {sim_concs_used}...')
            for sim_col in sim_cols:
                # sim_col is a tuple (self._SIMULATED, specie)
                sim_specie = sim_col[1]
                df[sim_col] = self._get_gaussian(x_range=x_range,
                                                 specie=sim_specie,
                                                 conc=sim_concs[sim_specie])
            df[self._SIM_ENV] = df[sim_cols].sum(axis=1)
            envelope = df[self._SIM_ENV]

        if experimental:
            _echo.echo(f'Processing experimental data...')
            df[self._RAW_EXP] = exp_data
            if self._CLEAN_EXP in df:
                df[self._CLEAN_EXP] = exp_data.copy()

        if gas_phase:
            _echo.echo(f'Calculating gas phase...')
            gas_gauss = self._get_gas_phase(x_range=x_range,
                                            exp_data=exp_data,
                                            gas_interval=gas_interval)
            df[self._GAS_PHASE] = gas_gauss
            _echo.echo(f'Removing gas phase from experimental...')
            df[self._CLEAN_EXP] -= df[self._GAS_PHASE]

        # Must go after simulate
        if autoscale:
            scale_factor = self._get_autoscale(sim_data=envelope,
                                               exp_data=exp_data)
            _echo.echo(f'Generating auto-scale factor {scale_factor}...')

        # Must go after autoscale
        if contaminate or decontaminate:
            contam_concs = {}
            for cont_col in contam_cols:
                # cont_col is a tuple (self._CONTAMINANT, specie)
                cont_specie = cont_col[1]
                contamination, contam_conc = self._get_contamination(
                    x_range=x_range,
                    contam_specie=cont_specie,
                    species=species,
                    scale_factor=scale_factor,
                    spectrum=contam_spectra[cont_specie]
                )
                contam_concs[cont_specie] = contam_conc
                df[cont_col] = contamination
                _echo.echo(f'Calculated contaminant concentrations at '
                           f'{contam_concs}...')

        if contaminate:
            _echo.echo('Adding contaminants to simulated data...')
            for cont_col in contam_cols:
                df[self._SIM_ENV] += df[cont_col]

        # Scale everything which needs to be scaled
        _echo.echo(f'Scaling data by factor {scale_factor}...')
        for sim_col in sim_cols:
            df[sim_col] *= scale_factor
        for cont_col in contam_cols:
            df[cont_col] *= scale_factor
        if self._SIM_ENV in df:
            df[self._SIM_ENV] *= scale_factor

        # Must happen after scaling
        if decontaminate:
            _echo.echo('Removing contaminants from experimental data...')
            for cont_col in contam_cols:
                df[self._CLEAN_EXP] -= df[cont_col]

        # Must go after scaling, gas_phase, and decontaminate
        if deconvolute:
            # Get a dict {specie: conc}, like sim_concs
            deconv_concs = self._get_deconvolution(
                x_range=x_range,
                species=deconv_species,
                species_guesses=sim_concs,
                scale_factor=scale_factor,
                to_fit=df[self._CLEAN_EXP]
            )
            for dec_col in deconv_cols:
                # dec_col is a tuple (self._DECONV, specie)
                dec_specie = dec_col[1]
                df[dec_col] = self._get_gaussian(x_range=x_range,
                                                 specie=dec_specie,
                                                 conc=deconv_concs[dec_specie])
            df[self._DECONV_ENV] = df[deconv_cols].sum(axis=1)
            _echo.echo(f'Deconvoluted experimental to get concentrations '
                       f'{deconv_concs}...')

            # More scaling
            for dec_col in deconv_cols:
                df[dec_col] *= scale_factor
            if self._DECONV_ENV in df:
                df[self._DECONV_ENV] *= scale_factor

        return df, scale_factor, deconv_concs, contam_concs

    def _get_autoscale(self,
                       sim_data: Union[pd.Series, None],
                       exp_data: Union[pd.Series, None]) -> float:
        """Gets the factor by which to automatically scale.

        Currently, gives 1.0 unless there's an experimental and a simulated,
        in which case it makes the peak experimental equal the peak envelope.

        Args:
            df: The pd.DataFrame with the gaussians to scale.
            sim_data: The simulated data you're going to scale.
            exp_data: Experimental data, if any, to autoscale to.
        """
        scale_factor = 1.0
        if exp_data is not None and sim_data is not None:
            scale_factor = max(exp_data) / max(sim_data)

        return scale_factor

    def _get_contamination(self, x_range: np.ndarray,
                           contam_specie: sym.Symbol,
                           species: List[sym.Symbol],
                           scale_factor: float,
                           spectrum: pd.Series, ) -> Tuple[np.ndarray, float]:
        """Get the contamination spectrum of the given specie.

        Args:
            x_range: The range you care about decontaminating.
            contam_specie: The specie you want the spetrum/concentration of.
            species: All of the species you want to be involved in the
                deconvolution of the contaminant spectrum. If there are many
                contaminants, put them all here.
            scale_factor: The scale factor of the simulated data. Keeps
                simulated and deconvoluted concentrations consistent. A
                placeholder for an intensity function.
            spectrum: The contaminant spectrum.

        Returns:
            The spectrum (a np.ndarray) and the concentration of the
            contaminant (a float).
        """

        # Deconvolute the spectra to get the conc of this contaminant
        deconved = self._get_deconvolution(
            x_range=np.asarray(spectrum.index),
            species=species,
            species_guesses={},
            scale_factor=scale_factor,
            to_fit=spectrum)
        contam_conc = deconved[contam_specie]

        # Make the gaussian of this contaminant on the x_range we care about.
        contam_gauss = self._get_gaussian(x_range=x_range,
                                          specie=contam_specie,
                                          conc=contam_conc)
        # TODO: I would put a ratio here, but I don't understand how it should
        # TODO: work, so I'll leave it out for now.

        return contam_gauss, contam_conc

    def _get_deconvolution(self, x_range: np.ndarray,
                           species: List[sym.Symbol],
                           species_guesses: Dict[sym.Symbol, float],
                           scale_factor: float,
                           to_fit: pd.Series) -> Dict[sym.Symbol, float]:
        """Get the deconvolution of experimental into species.

        Args:
            x_range: np.ndarray, the x-range on which to calculate gaussians.
            species: A list of sym.Symbols, the species for which we're 
                guessing concentrations.
            species_guesses: A dict {species: conc} of initial guesses.
                This should be sim_concs, as far as I can tell. It need not
                be populated.
            scale_factor: float, the factor by which to scale the gaussians.
                This is so that the concentrations of the simulated and
                deconvoluted data are comparable.
            to_fit: pd.Series, the spectrum to deconvolute.

        Returns:
            A dict of {species: conc} of the guessed concentration.
        """
        # Guess 1 for everything, unless there's simulated data, then use that.
        conc_guesses = np.ones(len(species))
        for index, specie in enumerate(species):
            if specie in species_guesses:
                conc_guesses[index] = species_guesses[specie]
            # In case we'd get an out of bounds error.
            conc_guesses[index] = max(0.0, conc_guesses[index])

        # Make the function to pass scipy.optimize.curve_fit. Needs signature
        # f(x-values, *y-values) -> ndarray of size to_fit.size
        def envelope_from_concs(x_range, *concs: float):
            envelope = np.zeros(x_range.size)
            for index, specie in enumerate(species):
                hill = self._get_gaussian(x_range, specie, concs[index])
                envelope += hill * scale_factor
            return envelope

        fitted_concs, covariance = optimize.curve_fit(f=envelope_from_concs,
                                                      xdata=x_range,
                                                      ydata=to_fit,
                                                      p0=conc_guesses,
                                                      bounds=(0, np.inf))
        return dict(zip(species, fitted_concs))

    def _get_gas_phase(self, x_range: np.ndarray,
                       exp_data: pd.Series,
                       gas_interval: Tuple[float, float]) -> np.ndarray:
        """Calculates the gas phase given an experimental.

        It picks as the peak of the gas phase the maximum of exp_data in
        gas_interval.

        Args:
            x_range: np.ndarray, the x-values on which to caluclate the
                gaussian.
            exp_data: pd.Series, the experimental data to use.
            gas_interval: An interval (2-ple), wherein you expect the peak of
                the gas phase to lie. Make it small, as otherwise it might pick
                the foot of a taller gaussian as the 'peak'.

        Returns:
            An np.ndarray of the gas phase gaussian.
        """
        # Get the range to search for the gas peak in.
        search = exp_data[gas_interval[0]:gas_interval[1]]

        # Get the location of the highest part of the experimental data in range
        peak_x = max(search.index, key=lambda index: search[index])
        peak = exp_data[peak_x]

        # Make a gaussian the same height as the experimental gas phase peak
        gas_gaussian = stats.norm.pdf(x_range, peak_x, self._SIGMA)
        gas_gaussian *= (peak / max(gas_gaussian))

        return gas_gaussian

    def _get_gaussian(self, x_range: np.ndarray,
                      specie: sym.Symbol,
                      conc: float) -> np.ndarray:
        """Calculates the (unscaled) gaussian for the given species.

        Args:
            x_range: np.ndarray, the x-values on which to caluclate the
                gaussian.
            specie: sym.Symbol, the specie for which to get the gaussian.
            conc: float, the concentration of the specie.

        Returns:
            A x_range-sized ndarray with a normal distribtuion centered at the
            binding energy of specie, with st.dev _SIGMA.
        """
        gaussian = np.zeros(x_range.size)

        for orbital in self.species_manager[specie].orbitals:
            hill = stats.norm.pdf(x_range, orbital.binding_energy, self._SIGMA)
            gaussian += conc * orbital.splitting * hill

        return gaussian

    def _get_x_range(self, species: List[sym.Symbol]) -> np.ndarray:
        """Picks an x-range on which to calculate gaussians.

        Args:
            species: A list of sym.Symbols, the species you want to be visible
                in the x_range.

        Returns:
            An np.ndarray on which we're going to calculate gaussians.
        """
        binding_energies = []
        for specie in species:
            for orbital in self.species_manager[specie].orbitals:
                binding_energies.append(orbital.binding_energy)

        x_lower = min(binding_energies) - self._PLOT_MARGIN
        x_upper = max(binding_energies) + self._PLOT_MARGIN
        x_range = np.arange(x_lower, x_upper, self._PLOT_RESOLUTION)

        return x_range

    # --- Plotting -----------------------------------------------------------

    def _plot(self, ax: plt.Axes,
              species: List[sym.Symbol],
              # Args passed to resample
              simulate: Optional[bool] = None,
              autoscale: Optional[bool] = None,  # Resample only
              experimental: Optional[bool] = None,
              gas_phase: Optional[bool] = None,
              decontaminate: Optional[bool] = None,  # TODO
              contaminate: Optional[bool] = None,  # TODO
              deconvolute: Optional[bool] = None,
              # Args passed only to plotter
              sim_gaussians: bool = True,
              envelope: bool = True,
              experimental_raw: bool = True,
              experimental_clean: bool = True,
              deconv_gaussians: bool = True,
              deconv_envelope: bool = True, ) -> plt.Axes:
        """Plot the XPS observable.

        Args:
            # TODO

        Returns:
            The plt.Axes of the plot.
        """
        # This makes an XPSObservable.
        xps_obs = self.resample(
            overwrite=False,
            species=species,
            simulate=simulate,
            autoscale=autoscale,
            experimental=experimental,
            gas_phase=gas_phase,
            decontaminate=decontaminate,
            contaminate=contaminate,
            deconvolute=deconvolute,
        )

        # Default behavior; yes I know there are shorter ways to do this, but
        # this is the most explicit.
        if simulate is None:
            simulate = True
        if experimental is None:
            experimental = True
        if gas_phase is None:
            gas_phase = True
        if deconvolute is None:
            deconvolute = True
        if decontaminate is None and contaminate is None:
            contaminants = True
        else:
            contaminants = decontaminate or contaminate

        # XPSExperiment.plot is overridden by Experiment.plot, but
        # this is an XPSObservable, so this doesn't make an infinite loop.
        return xps_obs.plot(ax=ax,
                            simulated=simulate,
                            sim_gaussians=sim_gaussians,
                            sim_envelope=envelope,
                            experimental=experimental,
                            exp_raw=experimental_raw,
                            exp_clean=experimental_clean,
                            gas_phase=gas_phase,
                            deconvoluted=deconvolute,
                            deconv_gaussians=deconv_gaussians,
                            deconv_envelope=deconv_envelope,
                            contaminants=contaminants)

    # --- Utility -------------------------------------------------------------

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict()
        d['species_concs'] = {str(symbol): conc for symbol, conc in
                              self._sim_concs.items()}
        d['species_manager'] = self.species_manager.as_dict()
        d['autoresample'] = self.autoresample
        d['autoscale'] = self.autoscale
        d['experimental'] = self._exp_data.to_json() if \
            self._exp_data is not None else None
        d['gas_interval'] = self._gas_interval
        d['scale_factor'] = self._scale_factor
        d['title'] = self.title
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        decode = monty.json.MontyDecoder().process_decoded
        d['species_concs'] = {sym.Symbol(name): conc for name, conc in
                              d['species_concs'].items()}
        d['species_manager'] = decode(d['species_manager'])
        if d['experimental'] is not None:
            d['experimental'] = pd.read_json(d['experimental'],
                                             typ='series',
                                             convert_axes=False)
            d['experimental'].index = d['experimental'].index.map(float)
        return cls(**d)


def simulate_xps(rsys: reaction.RxnSystem,
                 time: Optional[float] = None,
                 title: str = '',
                 species: List[sym.Symbol] = None,
                 ignore: List[sym.Symbol] = None,
                 x_range: Optional[np.ndarray] = None,
                 scale_factor: float = None,
                 experimental: pd.Series = None,
                 gas_interval: Tuple[float, float] = None,
                 contam_spectra: Optional[Dict[sym.Symbol, pd.Series]] = None,
                 deconv_species: Optional[List[sym.Symbol]] = None,
                 autoresample: bool = True,
                 autoscale: bool = True,
                 **options) -> XPSExperiment:
    """Simulate the given reaction system over time.

    Args:
        TODO: add rest
        rsys: ReactionsSystem, the reaction system to simulate
        time: The time until which to simulate.
        species: The species to include in the XPS.
        ignore: The species to not include in the XPS.
        autoresample: Decides if the XPS resamples on edits.
        autoscale: Decides if the XPS will automatically scale
            the gaussians and envelope to match the experimental data.
        experimental: The experimental value of the XPS.
        gas_interval: The interval in which the peak of the gas phase is
            in the XPS.
        scale_factor: The scale factor by which to scale the simulated
            gaussians in the XPS.
        title: The name to give the XPS, used in plotting.
        **options: Forwarded to scipy.integrate.solve_ivp

    Returns:
        A Solution object describing the solution.
    """
    if time is None:
        sol_t, sol_y = bulk_crn.solve_rsys_ode_till_eq(rsys, **options)
    else:
        sol_t, sol_y = bulk_crn.solve_rsys_ode(rsys, time, **options)

    cts = time_series.CRNTimeSeries(sol_t, sol_y, rsys)
    return cts.xps_with(title=title,
                        species=species,
                        ignore=ignore,
                        x_range=x_range,
                        scale_factor=scale_factor,
                        experimental=experimental,
                        gas_interval=gas_interval,
                        contam_spectra=contam_spectra,
                        deconv_species=deconv_species,
                        autoresample=autoresample,
                        autoscale=autoscale, )
