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

from __future__ import annotations
from typing import Dict, List, Optional, Tuple, Union

import copy
import json
import warnings

from matplotlib import pyplot as plt
import monty.json
import numpy as np
import pandas as pd
from scipy import integrate
from scipy import optimize
from scipy import stats
import sympy as sym

import dtcs
from dtcs.io import xps as xps_io
from dtcs.spec.model_input.relations import TPRateRelation
from dtcs.twin import twin_abc
from dtcs.spec import species
from dtcs.common import util
from dtcs import _logger
CRNTimeSeries: TypeAlias = 'CRNTimeSeries'

_PLOT_MARGIN = 5
_PLOT_RESOLUTION = 0.001
_SIGMA = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)
_NOISE_MAX = 0.005 # A proportion of the maximum value that will form the upper bound for added noise.

_REQUIRED = 'required'
_OPTIONAL = 'optional'


class XPSObservable(monty.json.MSONable):
    """TODO"""

    # Settings for the column names
    _SIMULATED = 'simulated'
    _EXPERIMENTAL = 'experimental'
    _CONTAMINANTS = 'contaminants'
    _DECONVOLUTED = 'deconvoluted'
    _SIM_ENV = (_SIMULATED, 'envelope')
    _EXP_CLEAN = (_EXPERIMENTAL, 'clean')
    _EXP_RAW = (_EXPERIMENTAL, 'raw')
    _GAS_PHASE = (_EXPERIMENTAL, 'gas_phase')
    _DECONV_ENV = (_DECONVOLUTED, 'envelope')

    _RESERVED = [_SIMULATED, _EXPERIMENTAL, _CONTAMINANTS, _DECONV_ENV,
                 _SIM_ENV[1], _EXP_CLEAN[1], _EXP_RAW[1], _GAS_PHASE[1],
                 _DECONV_ENV[1]]

    def __init__(self, df: pd.DataFrame,
                 species_manager: species.SpeciesManager,
                 title: str = ''):
        self.df: pd.DataFrame = df
        self.species_manager: species.SpeciesManager = species_manager
        self.title: str = title

    # --- Accessors ----------------------------------------------------------
    @property
    def x_range(self) -> Optional[np.ndarray]:
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
    def sim_envelope(self) -> Optional[pd.Series]:
        """The simulated envelope, the sum of all species' gaussians."""
        if self._SIM_ENV in self.df:
            return self.df[self._SIM_ENV]
        return None

    @property
    def sim_gaussians(self) -> Optional[pd.DataFrame]:
        """The simulated gaussians of the XPS observable."""
        if self.simulated is None:
            return None
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
    def experimental(self) -> Optional[pd.DataFrame]:
        """The experimental data. Might be None."""
        if self._EXPERIMENTAL in self.df:
            return self.df[self._EXPERIMENTAL]
        return None

    @experimental.setter
    def experimental(self, value):
        self.df[self._EXPERIMENTAL] = value
        return None

    @property
    def exp_clean(self) -> Optional[pd.Series]:
        """The experimental data without the gas phase or contaminants."""
        if self._EXP_CLEAN in self.df:
            return self.df[self._EXP_CLEAN]
        return None

    @property
    def exp_raw(self) -> Optional[pd.Series]:
        """The raw experimental data. Might be None."""
        if self._EXP_RAW in self.df:
            return self.df[self._EXP_RAW]
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
    def deconv_envelope(self) -> Optional[pd.Series]:
        """The envelope of the deconvoluted experimental data."""
        if self._DECONV_ENV in self.df:
            return self.df[self._DECONV_ENV]
        return None

    @property
    def deconv_gaussians(self) -> Optional[pd.DataFrame]:
        """The deconvoluted gaussians of the XPS observable."""
        if self.deconvoluted is None:
            return None
        # Gaussians have a species as their column name
        if self.deconvoluted is None:
            return None

        gauss_cols = []
        for col in self.deconvoluted:
            if isinstance(col, sym.Symbol):
                gauss_cols.append(col)
        if gauss_cols:
            return self.deconvoluted[gauss_cols]
        else:
            return None

    # --- Serialization ------------------------------------------------------
    def as_dict(self) -> dict:
        d = {
            '@module': self.__class__.__module__,
            '@class': self.__class__.__name__,
            '@version': dtcs.__version__,  # TODO: Better way to do this?
            'df': None,
            'species_manager': self.species_manager.as_dict(),
            'title': self.title,
        }

        def sanitize_column(column):
            category = column[0]
            name = util.symbol_to_name(column[1])

            return f'{category}\\{name}'

        df_copy = self.df.copy(deep=False)
        df_copy.columns = [sanitize_column(col) for col in df_copy.columns]
        # Pandas.to_dict() is couldn't do orient=table, so I'm using this
        #  roundabout method, sorry!
        d['df'] = json.loads(df_copy.to_json(
            orient='table',
        ))

        return d

    @classmethod
    def from_dict(cls, d: dict):
        decode = monty.json.MontyDecoder().process_decoded

        d['species_manager'] = decode(d['species_manager'])
        sm = d['species_manager']

        def desanitize_column(column):
            category, name = column.split('\\')
            symbol = sym.Symbol(name)
            if not (name in XPSObservable._RESERVED) and symbol in sm:
                name = symbol
            return category, name

        df = pd.read_json(
            json.dumps(d['df']),  # Frivolous string conversion
            orient='table',
            convert_axes=False,
        )
        columns = [desanitize_column(col) for col in df.columns]
        df.columns = pd.MultiIndex.from_tuples(columns)
        d['df'] = df

        return cls(**d)

    # --- Plotting -----------------------------------------------------------

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

        species = twin_abc._get_species_not_ignored(
            species,
            ignore,
            self.species_manager.symbols,
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
                    color = self.species_manager[specie].color
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
                color = self.species_manager[specie].color
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
                color = self.species_manager[specie].color
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
                color = self.species_manager[specie].color
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
        return np.sqrt(((self.exp_raw - self.sim_envelope) ** 2).mean())

    def mae(self):
        return np.sqrt((np.abs(self.exp_raw - self.sim_envelope)).mean())

    def correl(self):
        # TODO(Andrew): Move the following logic into a function somewhere,
        #  so you can access the 'best' experimental data with one line.
        if self.exp_clean is not None:
            exp_envelope = self.exp_clean
        elif self.exp_raw is not None:
            exp_envelope = self.exp_raw
        else:
            raise ValueError('Experimental data required!')

        return np.corrcoef(exp_envelope,
                           self.sim_envelope)[0][1]

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


class XPSExperiment(twin_abc.Experiment, XPSObservable):
    """A container for a simulated observable of an XPS experiment.

    Attributes:
        autoresample: Defaults to true, decides if it resamples on edits.
        autoscale: Defaults to true, decides if it will automatically scale the gaussians and envelope to match the experimental data.
        sim_concs: The concentrations of each species, for the creation of simulated data.
        species_manager: The SpeciesManager in use.
        ignore: Species to ignore during processing.
    """

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
                             f'sim_concs or experimental defined.')

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
                _logger.echo(f'No species specified to '
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
        self._scale_factor: float = scale_factor
        self._sim_concs: Dict[sym.Symbol, float] = sim_concs
        self._exp_input: Union[pd.Series, None] = experimental
        self._gas_interval: Union[Tuple[float, float], None] = gas_interval
        self._contam_concs: Dict[sym.Symbol, float] = None
        self._contam_spectra: Dict[sym.Symbol, pd.Series] = contam_spectra
        self._deconv_species: List[sym.Symbol] = deconv_species
        self._deconv_concs: Union[Dict[sym.Symbol, float], None] = None

        # Resample if autoresample is True.
        self._autoresample()

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

    # --- Mutators -----------------------------------------------------------
    # If autoresample is on, these will prompt an overwriting resample.

    def include(self, *species: sym.Symbol):
        """Add the species to the experiment, if they aren't there already.

        It won't simulate data for them (unless you also edit sim_concs),
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
            _logger.echo('Auto-resampling data...')
            self.resample(overwrite=True)

    # --- Calculations -------------------------------------------------------
    # Note that functions in this section will only modify state if
    # overwrite=True.

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
            gas_phase: Optional[bool] = False,
            gas_interval: Optional[Tuple[float, float]] = None,
            decontaminate: Optional[bool] = None,
            contaminate: Optional[bool] = None,
            contam_spectra: Optional[Dict[sym.Symbol, pd.Series]] = None,
            deconvolute: Optional[bool] = None,
            deconv_species: Optional[List[sym.Symbol]] = None,
            augment = False,
    ) -> XPSObservable:
        """TODO"""
        # --- Parse Arguments ------------------------------------------------
        # Handle: species, ignore
        species = twin_abc._get_species_not_ignored(species, ignore,
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
            exp_data = self._exp_input  # Might still be None

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
                deconvolute = _logger.prompt_yn('Deconvolute experimental?',
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
        # x_range might still be None, in which case _resample() will scale it
        # automatically.

        # --- Resample and Overwrite -----------------------------------------
        df, scale_factor, deconv_concs, contam_concs = self._do_resample(
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
            augment=augment,
        )

        if overwrite:
            self._scale_factor = scale_factor
            self._sim_concs = sim_concs
            self._exp_input = exp_data
            self._gas_interval = gas_interval
            self._contam_concs = contam_concs
            self._contam_spectra = contam_spectra
            self._deconv_concs = deconv_concs
            self.df = df
        return XPSObservable(species_manager=self.species_manager,
                             df=df,
                             title=self.title)

    def _do_resample(
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
            augment: bool,
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
        # TODO(Andrew): Fix these assertions
        # assert not ((decontaminate or contaminate) and not contam_spectra)
        # assert not ((decontaminate or contaminate) and not contam_species)
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
            x_range = XPSExperiment._get_x_range(
                species_manager=self.species_manager,
                species=species_to_plot,
            )
            _logger.echo(f'Using automatically-generated x-range '
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
            columns.append(self._EXP_RAW)
            # Gas phase and decontaminate will clean up the raw experimental.
        if gas_interval or decontaminate:
            columns.append(self._EXP_CLEAN)
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

        df = XPSExperiment._get_dataframe(x_range=x_range,
                                          columns=columns)

        # --- Calculations ---------------------------------------------------
        if simulate:
            if _logger.do_echo:
                sim_concs_used = {specie: conc for specie, conc in
                                  sim_concs.items() if specie in sim_species}
                _logger.echo(f'Simulating XPS experiment for species with '
                           f'concentrations {sim_concs_used}...')
            for sim_col in sim_cols:
                # sim_col is a tuple (self._SIMULATED, specie)
                sim_specie = sim_col[1]
                df[sim_col] = XPSExperiment._get_gaussian(
                    species_manager=self.species_manager,
                    x_range=x_range,
                    specie=sim_specie,
                    conc=sim_concs[sim_specie],
                )
            df[self._SIM_ENV] = df[sim_cols].sum(axis=1)
            envelope = df[self._SIM_ENV]

        if experimental:
            _logger.echo(f'Processing experimental data...')
            df[self._EXP_RAW] = exp_data
            if self._EXP_CLEAN in df:
                df[self._EXP_CLEAN] = exp_data.copy()

        if gas_phase:
            _logger.echo(f'Calculating gas phase...')
            gas_gauss = self._get_gas_phase(x_range=x_range,
                                            exp_data=exp_data,
                                            gas_interval=gas_interval)
            df[self._GAS_PHASE] = gas_gauss
            _logger.echo(f'Removing gas phase from experimental...')
            df[self._EXP_CLEAN] -= df[self._GAS_PHASE]

        # Must go after simulate
        if autoscale:
            # Andrew: This line here was used for intelligent autoscaling.
            #  species_info = [self.species_manager[specie] for specie in species]
            scale_factor = self._get_autoscale(sim_envelope=envelope,
                                               exp_envelope=exp_data)
            _logger.echo(f'Generating auto-scale factor {scale_factor}...')

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
                _logger.echo(f'Calculated contaminant concentrations at '
                           f'{contam_concs}...')

        if contaminate:
            _logger.echo('Adding contaminants to simulated data...')
            for cont_col in contam_cols:
                df[self._SIM_ENV] += df[cont_col]

        # Scale everything which needs to be scaled
        _logger.echo(f'Scaling data by factor {scale_factor}...')
        for sim_col in sim_cols:
            df[sim_col] *= scale_factor
        for cont_col in contam_cols:
            df[cont_col] *= scale_factor
        if self._SIM_ENV in df:
            df[self._SIM_ENV] *= scale_factor

        # Must happen after scaling
        if decontaminate:
            _logger.echo('Removing contaminants from experimental data...')
            for cont_col in contam_cols:
                df[self._EXP_CLEAN] -= df[cont_col]

        # Must go after scaling, gas_phase, and decontaminate
        if deconvolute:
            # Get a dict {specie: conc}, like sim_concs
            deconv_concs = self._get_deconvolution(
                x_range=x_range,
                species=deconv_species,
                species_guesses=sim_concs,
                scale_factor=scale_factor,
                to_fit=df[self._EXP_CLEAN]
            )
            for dec_col in deconv_cols:
                # dec_col is a tuple (self._DECONV, specie)
                dec_specie = dec_col[1]
                df[dec_col] = XPSExperiment._get_gaussian(
                    species_manager=self.species_manager,
                    x_range=x_range,
                    specie=dec_specie,
                    conc=deconv_concs[dec_specie]
                )
            df[self._DECONV_ENV] = df[deconv_cols].sum(axis=1)
            _logger.echo(f'Deconvoluted experimental to get concentrations '
                       f'{deconv_concs}...')

            # More scaling
            for dec_col in deconv_cols:
                df[dec_col] *= scale_factor
            if self._DECONV_ENV in df:
                df[self._DECONV_ENV] *= scale_factor

        if augment:
            m = max(df[self._SIM_ENV])
            df[self._SIM_ENV] += (_NOISE_MAX*m)*np.random.normal(0,1,len(df[self._SIM_ENV]))

        return df, scale_factor, deconv_concs, contam_concs

    # --- Calculation Helper-functions ---------------------------------------

    @staticmethod
    def _get_dataframe(x_range: np.ndarray,
                       columns: List) -> pd.DataFrame:
        """Get a properly-formatted dataframe with the given x_range and
        columns."""
        row_index = pd.Index(x_range, name='eV')
        col_index = pd.MultiIndex.from_tuples(columns)
        return pd.DataFrame(data=0, index=row_index, columns=col_index)

    @staticmethod
    def _get_autoscale(sim_envelope: Union[pd.Series, None],
                       exp_envelope: Union[pd.Series, None]) -> float:
        """Gets the factor by which to automatically scale.

        Currently, gives 1.0 unless there's an experimental and a simulated,
        in which case it makes the two envelopes match at the location of the
        highest peak.

        Args:
            df: The pd.DataFrame with the gaussians to scale.
            sim_envelope: The simulated data you're going to scale.
            exp_envelope: Experimental data, if any, to autoscale to.
            species: The species for which the highest peak will be made to
                match.
        """
        scale_factor = 1.0
        if exp_envelope is not None and sim_envelope is not None:
            assert exp_envelope.index.equals(sim_envelope.index), 'Simulated' \
                    'and experimental indices should be the same by now.'

            lsq_sol, residuals, _, _ = np.linalg.lstsq(
                a=np.asarray([sim_envelope]).T,  # coefficient matrix, A
                b=exp_envelope,  # goal vector, b
            )  # This should do argmin_x(|Ax-b|^2), which is 1-dimensional least squares (x is a 1-by-1 matrix)

            scale_factor = lsq_sol[0]
            # TODO(Andrew) Maybe add this as a debugging utility, that would
            #  be cool.
            # dtcs._debug_exports = {
            #     'peaks': peaks,
            #     'exp_peaks': exp_envelope.loc[peaks],
            #     'sim_peaks': sim_envelope.loc[peaks],
            # }
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
        contam_gauss = XPSExperiment._get_gaussian(
            species_manager=self.species_manager,
            x_range=x_range,
            specie=contam_specie,
            conc=contam_conc
        )
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
                hill = self._get_gaussian(species_manager=self.species_manager,
                                          x_range=x_range,
                                          specie=specie,
                                          conc=concs[index])
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
        gas_gaussian = stats.norm.pdf(x_range, peak_x, _SIGMA)
        gas_gaussian *= (peak / max(gas_gaussian))

        return gas_gaussian

    @staticmethod
    def _get_x_range(species_manager: species.SpeciesManager,
                     species: List[sym.Symbol]) -> np.ndarray:
        """Picks an x-range on which to calculate gaussians.

        Args:
            species: A list of sym.Symbols, the species you want to be visible
                in the x_range.

        Returns:
            An np.ndarray on which we're going to calculate gaussians.
        """
        binding_energies = []
        for specie in species:
            for orbital in species_manager[specie].orbitals:
                binding_energies.append(orbital.binding_energy)

        x_lower = min(binding_energies or [0]) - _PLOT_MARGIN
        x_upper = max(binding_energies or [0]) + _PLOT_MARGIN
        x_range = np.arange(x_lower, x_upper, _PLOT_RESOLUTION)

        return x_range

    @staticmethod
    def _get_gaussian(species_manager: species.SpeciesManager,
                      x_range: np.ndarray,
                      specie: sym.Symbol,
                      conc: float) -> np.ndarray:
        """Calculates the (unscaled) gaussian for the given species.

        Args:
            species_manager: The species manager to use to get the binding
                energy.
            x_range: np.ndarray, the x-values on which to caluclate the
                gaussian.
            specie: sym.Symbol, the specie for which to get the gaussian.
            conc: float, the concentration of the specie.

        Returns:
            A x_range-sized ndarray with a normal distribtuion centered at the
            binding energy of specie, with st.dev _SIGMA.
        """
        gaussian = np.zeros(x_range.size)

        for orbital in species_manager[specie].orbitals:
            hill = stats.norm.pdf(x_range, orbital.binding_energy, _SIGMA)
            gaussian += conc * orbital.splitting * hill

        return gaussian

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
              deconvolute: Optional[bool] = False,
              augment: bool = False,
              # Args passed only to plotter
              # sim_gaussians: bool = True,
              # envelope: bool = True,
              # experimental_raw: bool = True,
              # experimental_clean: bool = True,
              # deconv_gaussians: bool = True,
              # deconv_envelope: bool = True,
              **kwargs,
    ) -> plt.Axes:
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
            augment=augment,
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
                            # sim_gaussians=sim_gaussians,
                            # sim_envelope=envelope,
                            experimental=experimental,
                            # exp_raw=experimental_raw,
                            # exp_clean=experimental_clean,
                            gas_phase=gas_phase,
                            deconvoluted=deconvolute,
                            # deconv_gaussians=deconv_gaussians,
                            # deconv_envelope=deconv_envelope,
                            contaminants=contaminants,
                            **kwargs)

    # --- Utility -------------------------------------------------------------

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict()
        d['sim_concs'] = {str(symbol): conc for symbol, conc in
                              self._sim_concs.items()}
        d['species_manager'] = self.species_manager.as_dict()
        d['autoresample'] = self.autoresample
        d['autoscale'] = self.autoscale
        d['experimental'] = self._exp_input.to_json() if \
            self._exp_input is not None else None
        d['gas_interval'] = self._gas_interval
        d['scale_factor'] = self._scale_factor
        d['title'] = self.title
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        decode = monty.json.MontyDecoder().process_decoded
        d['sim_concs'] = {sym.Symbol(name): conc for name, conc in
                              d['sim_concs'].items()}
        d['species_manager'] = decode(d['species_manager'])
        if d['experimental'] is not None:
            d['experimental'] = pd.read_json(d['experimental'],
                                             typ='series',
                                             convert_axes=False)
            d['experimental'].index = d['experimental'].index.map(float)
        return cls(**d)


class XPSSolutionSystem:
    """Process and visualize a system of experiment

    This system performs automatic global scaling as well as other processing.

    Instance Variables:
        ignore: A list of sympy symbols to ignore
        systems: A list of lists of XPS experiments.
        multipliers: A list of multipliers used in simulation.
    """
    def __init__(self, experiments: List, time_series: List, experimental_files: List[str], multipliers: List[float]):
        """Create a new solution system

        Given a list of lists of XPSExperiments (for all variations that are simulated) and a
        corresponding list of experimental files, experimental data is added to the XPSExperiment
        objects and the data is then auto-scaled.
        """
        assert len(experiments) > 0, 'Must pass at least one experiment to a Solution System.'
        assert len(time_series) == len(experiments), 'Must pass the same number of simulations and time series.'

        self.systems = experiments
        self.time_series = time_series

        # Iterate through XPS classes, adding experimental data.
        for i, f in enumerate(experimental_files):
            xps_exp = xps_io.read_new_data(f)[0]
            series = pd.Series(data=xps_exp.intensities, index=xps_exp.binding_energies)
            for s in self.systems[i]:
                s.experimental = series

        scaling_factor, max_index = self.max_scaling_factor()
        self.scale(scaling_factor)
        self._ignore = []
        self.multipliers = multipliers

        print('scaling factor:', scaling_factor, '\tmax index:', max_index)

    def _process(self):
        """Rescale solution data.
        """
        scaling_factor, max_index = self.max_scaling_factor()
        self.scale(scaling_factor)

    @property
    def ignore(self) -> List[sym.Symbol]:
        return self._ignore

    @ignore.setter
    def ignore(self, ignore: List[sym.Symbol]) -> None:
        self._ignore = ignore
        for sys in self.systems:
            for s in sys:
                s.ignore = self._ignore
        self._process()

    def max_scaling_factor(self) -> int:
        """Find the maximum peak and return the scaling factor for that data
        and the system index."""
        # If there is no experimental data, do not rescale
        if self.systems[0][0].exp_input is None:
            return 1, 0

        self.systems[0][0].scale_factor = 1
        max_env = max(self.systems[0][0].sim_envelope)
        max_index = -1
        max_exp = max(self.systems[0][0].exp_input)

        # TODO: I changed this as little as possible but it seems to put
        # TODO: envelope to be larger than peak.
        # Find max peak of systems
        for i in range(len(self.systems)):
            for j in range(len(self.systems[i])):
                self.systems[i][j].scale_factor = 1
                potential = max(self.systems[i][j].experimental)
                if potential > max_exp:
                    max_exp = potential
                    max_index = i
                    max_env = max(self.systems[i][j].sim_envelope)

        return max_exp / max_env, max_index

    def scale(self, scaling_factor):
        """Scale experimental data intensity to match that of the experimental data.
        """
        for sys in self.systems:
            for s in sys:
                s.scale_factor = scaling_factor

    def __getitem__(self, index):
        return self.systems[index]

    def time_series_at(self, sys: int, exp: int) -> CRNTimeSeries:
        return self.time_series[sys][exp]

    def plot(self, index):
        cols = len(self.multipliers)
        rows = len(self.systems[index]) // cols
        fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=(60, 60))
        plt.figure(fig.number)

        for j, ax in enumerate(fig.axes):
            plt.sca(ax)
            self.systems[index][j].plot()
            plt.title(f'Eq: {str(j // cols)} Const: {str(j % cols)}')
        plt.show()


class XPSInitializationData:
    """Stores the reaction constants and raw experimental file location.

    This class is simply syntactic sugar for a map.

    Instance Variables:
        name: The name of this specific experiment.
        constants: A list reaction constants for this experiment, with an order corresponding to that of the reaction system to be simulated.
        temp: The temperature associated with this experiment.
        pressure: The pressure associated with this experiment.
        experimental_file: The location of the experimental XPS file to be used for comparison
    """
    def __init__(self, title: str, temp: float, pressure: float, experimental_file: str = None, constants: List[float] = None):
        self.title = title
        self.constants = constants
        self.experimental_file = experimental_file
        self.temp = temp
        self.pressure = pressure


class XPSSystemRunner:
    """Orchestrate an end-to-end XPS analysis pipeline.

    Given a reaction system and initialization data objects, runs the specified simulations, returning
    a solution system object that can be used to visualize/analyze results.

    Instance Variables:
        rsys_generator: A function that takes in a list of reaction constants and returns a reaction
        system to be used.
        time: The amount of time for which the simulation is run.
        initializer_data: A list of XPSInitializationData objects that specify the various
        simulations to be run.
        multipliers: A list of reaction constant multipliers to be used.
        solutions: A list of lists of simulated data for an XPS experiment.
        complete: A list of boolean values representing whether all systems have been simulated.
        tp_correlator: An optional temperature pressure corrl
    """

    def __init__(self,
                 rsys_generator,
                 time: float,
                 initializer_data: List[XPSInitializationData],
                 multipliers: List[float],
                 tp_correlator: Optional[TPRateRelation] = None):
        """Create a new XPS system runner.
        """
        if len(multipliers) == 0:
            multipliers = [1]

        self.rsys_generator = rsys_generator
        self.time = time
        self.multipliers = multipliers
        self.solutions = [None for _ in range(len(initializer_data))]
        self.time_series = [None for _ in range(len(initializer_data))]
        self.complete = [False for _ in range(len(initializer_data))]
        self.tp_correlator = tp_correlator


        # Auto-scale reaction constants if applicable
        if tp_correlator is None:
            self.initializer_data = initializer_data
        else:
            self.generate_constants_from_correlation(initializer_data)

    def generate_constants_from_correlation(self, initializer_data: List[XPSInitializationData]):
        """Given a list of initialization data, and a tp_correlator, determine reaction constants.
        """
        tpc = self.tp_correlator
        exp: Optional[XPSInitializationData] = None
        # Find the reaction with which the initial data in the correlator corresponds
        for d in initializer_data:
            if d.pressure == tpc.init_pressure and d.temp == tpc.init_temp:
                exp = d

        if exp is None:
            raise ValueError("The TP correlators' initial temperature and pressure do not "
                             + "correspond with any experiment in the initialization data")

        for d in initializer_data:
            d.constants = tpc.constants(d.temp, d.pressure)
            print(d.constants)

        self.initializer_data = initializer_data


    def simulate(self, index: int):
        data = self.initializer_data[index]
        sols = []
        cts = []
        sol_mult_1 = None
        cts_mult_1 = None

        for i in range(len(data.constants)):
            for j in range(len(self.multipliers)):
                if self.multipliers[j] == 1 and sol_mult_1:
                    sols.append(sol_mult_1)
                    cts.append(cts_mult_1)

                else:
                    scaled = list(data.constants)
                    scaled[i] *= self.multipliers[j]

                    rsys = self.rsys_generator(scaled)
                    # TODO: simulate_xps_with_cts old?
                    s, ts = simulate.simulate_xps_with_cts(rsys, time=self.time, title=data.title + " Eq: "
                        + str(i) + "Constant: " + str(j))

                    cts.append(ts)
                    sols.append(s)
                    if self.multipliers[j] == 1:
                        sol_mult_1 = s
                        cts_mult_1 = ts

                print('Solved for ('+str(i)+', '+str(j)+')')
                print(scaled)
                print('\n')
        self.solutions[index] = sols
        self.time_series[index] = cts
        self.complete[index] = True

    def system(self) -> XPSSolutionSystem:
        """Returns a solution system if all experiments have been simulated.
        """
        if not all(self.complete):
            raise AssertionError("All experiments have not yet been simulated.")
        return XPSSolutionSystem(self.solutions, self.time_series, [x.experimental_file for x in self.initializer_data if x.experimental_file],
                                 self.multipliers)


def simulate_xps(*args, **kwargs):
    raise NotImplementedError()
