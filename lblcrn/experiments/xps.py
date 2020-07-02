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
from lblcrn import _echo


class XPSObservable:
    """A wrapper around a dataframe, formatted to act like an XPS observable.

    Attributes:
        df: The pandas DataFrame in which the observables are stored.
        title: A str used for the title during plotting.
    """

    _COLORS = ['red', 'green', 'orange', 'blue', 'purple', 'pink', 'yellow',
               'brown', 'cyan']

    # Settings for the column names
    _SIMULATED = 'simulated'
    _EXPERIMENTAL = 'experimental'
    _CONTAMINANT = 'contaminant'
    _ENVELOPE = (_SIMULATED, 'envelope')
    _GAS_PHASE = (_EXPERIMENTAL, 'gas_phase')
    _DECON_ENV = (_EXPERIMENTAL, 'decon_envelope')
    _CLEAN_EXP = (_EXPERIMENTAL, 'clean')
    _RAW_EXP = (_EXPERIMENTAL, 'raw')

    def __init__(self, species_manager: species.SpeciesManager,
                 df: pd.DataFrame, title: str = ''):
        self.species_manager = species_manager
        self.df = df
        self.title = title

    @property
    def envelope(self) -> Optional[pd.Series]:
        """The simulated envelope, the sum of all species' gaussians."""
        if self._ENVELOPE in self.df:
            return self.df[self._ENVELOPE]
        return None

    @property
    def experimental_raw(self) -> Optional[pd.Series]:
        """The raw experimental data. Might be None."""
        if self._RAW_EXP in self.df:
            return self.df[self._RAW_EXP]
        return None

    @property
    def experimental_clean(self) -> Optional[pd.Series]:
        """The experimental data without the gas phase or contaminants."""
        if self._CLEAN_EXP in self.df:
            return self.df[self._CLEAN_EXP]
        return None

    @property
    def decon_envelope(self):
        """The envelope of the deconvoluted experimental data."""
        if self._DECON_ENV in self.df:
            return self.df[self._DECON_ENV]
        return None

    @property
    def experimental(self) -> Optional[pd.DataFrame]:
        """The experimental data. Might be None."""
        if self._EXPERIMENTAL in self.df:
            return self.df[self._EXPERIMENTAL]
        return None

    @property
    def simulated(self) -> Optional[pd.DataFrame]:
        """The simulated data. Might be None."""
        if self._SIMULATED in self.df:
            return self.df[self._EXPERIMENTAL]
        return None

    @property
    def gas_phase(self) -> Optional[pd.Series]:
        """The gas phase part of the spectrum. Might be None."""
        if self._GAS_PHASE in self.df:
            return self.df[self._GAS_PHASE]
        return None

    @property
    def gaussians(self) -> Optional[pd.DataFrame]:
        """The gaussians of the XPS observable."""
        # Columns are a pair (simulated/experimental, str/species)
        # If they're a species, we want them.
        gauss_cols = []
        for col in self.df.columns:
            if isinstance(col[1], sym.Symbol):
                gauss_cols.append(col)
        if gauss_cols:
            return self.df[gauss_cols]
        else:
            return None

    @property
    def gaussians_deconvoluted(self) -> Optional[pd.DataFrame]:
        """The deconvoluted gaussians of the XPS observable."""
        # Columns are a pair (simulated/experimental, str/species)
        # If they're a (experimental, species), we want them.
        gauss_cols = []
        for col in self.df.columns:
            if col[0] == self._EXPERIMENTAL and isinstance(col[1], sym.Symbol):
                gauss_cols.append(col)
        if gauss_cols:
            return self.df[gauss_cols]
        else:
            return None

    @property
    def gaussians_simulated(self) -> Optional[pd.DataFrame]:
        """The simulated gaussians of the XPS observable."""
        # If they picked simulated gaussians, return those
        gauss_cols = []
        for col in self.df.columns:
            if col[0] == self._SIMULATED and isinstance(col[1], sym.Symbol):
                gauss_cols.append(col)
        if gauss_cols:
            return self.df[gauss_cols]
        else:
            return None

    @property
    def x_range(self) -> np.ndarray:
        """The x-values, energies, on which there is data."""
        return np.asarray(self.df.index)

    def plot(self, ax: plt.Axes = None,
             # TODO: perhaps rename these all to show_ or show=?
             simulate: bool = True,
             sim_gaussians: bool = True,
             envelope: bool = True,
             experimental: bool = True,
             experimental_raw: bool = True,
             experimental_clean: bool = True,
             gas_phase: bool = True,
             deconvolute: bool = True,
             decon_gaussians: bool = True,
             decon_envelope: bool = True,
             **kwargs) -> plt.Axes:
        """Default plotting behavior for an XPS observable.

        Args:
            ax: The plt.Axes on which to plot.
            **kwargs: Forwarded. # TODO: how to use?
            # TODO

        Returns:
            # TODO
        """
        if ax is None:
            ax = plt.gca()

        # Sort the gaussian columns so shorter ones show in front.
        def gauss_col_sort_key(col):
            return max(self.df[col])

        # Some switches disable more than one thing
        if not simulate:
            sim_gaussians = False
            envelope = False
        if not experimental:
            experimental_raw = False
            experimental_clean = False
            gas_phase = False
            deconvolute = False
        if not deconvolute:
            decon_gaussians = False
            decon_envelope = False

        if sim_gaussians and self.gaussians_simulated is not None:
            for index, column in enumerate(sorted(self.gaussians_simulated,
                                                  key=gauss_col_sort_key,
                                                  reverse=True)):
                ax.fill(self.x_range, self.df[column],
                        label=f'Sim. {column[1]}',
                        color=self._COLORS[index])

        if decon_gaussians and self.gaussians_deconvoluted is not None:
            for index, column in enumerate(sorted(self.gaussians_deconvoluted,
                                                  key=gauss_col_sort_key,
                                                  reverse=True)):
                ax.plot(self.x_range, self.df[column],
                        label=f'Deconv. {column[1]}',
                        color=self._COLORS[index])

        if decon_envelope and self.decon_envelope is not None:
            ax.plot(self.x_range, self.decon_envelope,
                    label='Deconv. Envelope', color='black', linewidth=3)

        if gas_phase and self.gas_phase is not None:
            ax.fill(self.x_range, self.gas_phase, label='Gas Phase',
                    color='#666666')

        if envelope and self.envelope is not None:
            ax.plot(self.x_range, self.envelope, label='Sim. Envelope',
                    color='black', linewidth=3)

        if experimental_raw and self.experimental_raw is not None:
            ax.plot(self.x_range, self.experimental_raw,
                    label='Raw Experimental', color='#AAAAAA', linewidth=3)

        if experimental_clean and self.experimental_clean is not None:
            ax.plot(self.x_range, self.experimental_clean,
                    label='Clean Experimental', color='#DDAADD', linewidth=3)  # TODO: I need more colors/line styles

        ax.legend()
        ax.set_title(self.title)
        ax.invert_xaxis()  # XPS Plots are backwards
        return ax

    # --- Data Analysis ------------------------------------------------------

    def rmse(self):
        return np.sqrt(metrics.mean_squared_error(self.experimental_raw,
                                                  self.envelope))

    def mae(self):
        return metrics.mean_absolute_error(self.experimental_raw,
                                           self.envelope)

    def integral_diff_outside_experimental(self):
        return integrate.trapz(self.envelope - self.experimental_raw,
                               self.x_range)

    def integral_diff_inside_experimental(self):
        return integrate.trapz(self.experimental_raw - self.envelope,
                               self.x_range)

    def integral_diff_between(self):
        return integrate.trapz(np.absolute(self.experimental_raw -
                                           self.envelope), self.x_range)

    def envelope_integral(self):
        return integrate.trapz(self.envelope, self.x_range)


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
                 decon_species: Optional[List[sym.Symbol]] = None,
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
        if decon_species is None:
            decon_species = []
        if experimental is not None:
            experimental = experimental.copy()
        if species is None:
            species = []
            species.extend(sim_concs.keys())
            species.extend(contam_spectra.keys())
            species.extend(decon_species)
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
        self._decon_species: List[sym.Symbol] = decon_species

        # Resample if autoresample is True.
        self._autoresample()

    # --- Accessors ----------------------------------------------------------

    @property
    def gas_interval(self) -> Optional[Tuple[float, float]]:
        """The interval in which the gas phase peak should be. May be None."""
        return self._gas_interval

    @property
    def scale_factor(self) -> float:
        """The factor by which to scale the simulated data."""
        return self._scale_factor

    @property
    def species(self) -> List[sym.Symbol]:
        """All the species in the experiment."""
        return self._species

    # --- Mutators -----------------------------------------------------------
    # If autoresample is on, these will prompt an overwriting resample.

    def set_experimental(self, experimental: pd.Series):
        """Set the experimental data, prompting an autoresample."""
        self._exp_data = experimental
        self._autoresample()

    def del_experimental(self):
        """Deletes experimental, prompting an autoresample.."""
        self._exp_data = None
        del self.df['experimental']
        del self.df['gas_phase']
        self._autoresample()

    def set_gas_interval(self, lower: float, upper: float):
        """Sets the interval in which the gas phase peak should be.

        Give bounds for where the _peak_ is: if your bounds are too broad, it
        will get greedy, as it assumes that that interval is dominated by the
        gas phase. You can make lower equal to upper if you know the peak.

        Args:
            lower: float, the lower bound.
            upper: float, the upper bound.
        """
        if lower > upper:
            raise ValueError(f'Invalid interval ({lower}, {upper})')
        self._gas_interval = (lower, upper)
        self._autoresample()

    def del_gas_interval(self):
        """Deletes the gas interval, prompting a resample."""
        self._gas_interval = None
        del self.df['gas_phase']
        self._autoresample()

    def scale(self, factor: Optional[float] = None,
              by: Optional[float] = None,
              to: Optional[float] = None) -> float:
        """Scales according to the inputs you give.

        If you give no inputs, it autoscales. If you give inputs, it disables
        autoscaling.

        Args:
            factor: Optional, float, sets the scale factor exactly.
            by: Optional, float, multiplies the scale factor by this amount.
            to: Optional, float, makes the peak simulated envelope equal to it.

        Returns:
            The new scale factor.

        Raises:
            ValueError: If you supply more than one argument.
        """
        if sum([factor is None, by is None, to is None]) > 1:
            raise ValueError('Only one or zero arguments is acceptable.')

        if factor is not None:
            self.autoscale = False
            self._scale_factor = factor
        elif by is not None:
            self.autoscale = False
            self._scale_factor *= by
        elif to is not None:
            self.autoscale = False
            raise NotImplementedError()  # TODO
        else:
            raise NotImplementedError()  # TODO

        self._autoresample()
        return self._scale_factor

    def include(self, *species: sym.Symbol):
        """Add the species to the experiment.

        It won't simulate data for them (unless you also edit species_concs),
        but it will involve the species in deconvolution and the like.

        Prompts and autoresample.
        """
        self._species.extend(species)
        self._autoresample()

    def drop(self, *species: sym.Symbol):
        """Remove the species from the experiment.

        It won't deconvolute with or simulate data from them.

        Prompts and autoresample.
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
            decon_species: Optional[List[sym.Symbol]] = None,
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
            decon_species:

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
        species = self._get_species_not_ignored(species, ignore)

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
            raise ValueError('Decontamination requires experimental.')
        elif contaminate and not experimental:
            raise ValueError('Contamination requires experimental.')
        # (de)contaminate requires contam_species to be nonempty
        elif decontaminate and not contam_species:
            raise ValueError('Decontamination requries contaminants.')
        elif contaminate and not contam_species:
            raise ValueError('Contamination requires contaminants.')
        else:
            # decontaminate defaults to experimental, contaminate to False
            if decontaminate is None:
                decontaminate = experimental
            if contaminate is None:
                contaminate = False

        # Handle: decon_species
        if decon_species is None:
            decon_species = self._decon_species
            # If the default is unspecified, default to species
            if not decon_species:
                decon_species = species
        # Only take into account species given by _get_species_not_ignored
        decon_species = [specie for specie in decon_species
                         if specie in species]  # Everything subsets species

        # Handle: deconvolute
        # deconvolute requires experimental
        if deconvolute and not experimental:
            raise ValueError('Deconvolution requires experimental.')
        # deconvolute requires decon_species to be nonempty
        elif deconvolute and not decon_species:
            raise ValueError('Deconvolution requires more than zero species.')
        elif deconvolute is None:
            # Default deconvolute to experimental
            deconvolute = experimental

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
        df = self._resample(
            x_range=x_range,
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
            decon_species=decon_species,
        )

        if overwrite:
            self._x_range = x_range
            self._sim_concs = sim_concs
            self._scale_factor = scale_factor
            self._exp_data = exp_data
            self._gas_interval = gas_interval
            self._contam_spectra = contam_spectra
            self.df = df
        return XPSObservable(species_manager=self.species_manager,
                             df=df,
                             title=self.title)

    def _resample(
            self,
            # Input Data
            x_range: Union[np.ndarray, None],
            sim_concs: Dict[sym.Symbol, float],
            sim_species: List[sym.Symbol],
            scale_factor: float,
            exp_data: pd.Series,
            gas_interval: Tuple[float, float],
            contam_spectra: Dict[sym.Symbol, pd.Series],
            contam_species: List[sym.Symbol],
            decon_species: List[sym.Symbol],
            # Flags
            simulate: bool,
            autoscale: bool,
            experimental: bool,
            gas_phase: bool,
            decontaminate: bool,
            contaminate: bool,
            deconvolute: bool,
    ) -> pd.DataFrame:
        """TODO"""
        # Find the x_range based on what species we need to see.
        if x_range is None:
            species_to_plot = []
            if simulate:
                species_to_plot.extend(sim_species)
            # Don't bother with deconv_/contam_species, because if they're
            # defined, then the x_range will be exp_data.index anyway.
            x_range = self._get_x_range(species=species_to_plot)

        # Make the column (pd.MultiIndex) and row (pd.Index) indices.
        columns = []
        if simulate:
            columns.append(self._ENVELOPE)
            sim_cols = [(self._SIMULATED, specie) for specie in sim_species]
            columns.extend(sim_cols)
        if experimental:
            columns.append(self._RAW_EXP)
            # Gas phase and decontaminate will clean up the raw experimental.
            if gas_phase or decontaminate:
                columns.append(self._CLEAN_EXP)
            if gas_phase:
                columns.append(self._GAS_PHASE)
            if deconvolute:
                decon_cols = [(self._EXPERIMENTAL, specie)
                              for specie in decon_species]
                columns.extend(decon_cols)
            # TODO: (Andrew) Add a decontaminated spectrum
            # TODO: It would be the extrapolated contribution from contaminants.
            # TODO: We'd calculate it while decontaminating.

        # Build the DataFrame
        row_index = pd.Index(x_range, name='eV')
        col_index = pd.MultiIndex.from_tuples(columns)
        df = pd.DataFrame(data=0, index=row_index, columns=col_index)

        if experimental:
            df[self._RAW_EXP] = exp_data
            df[self._CLEAN_EXP] = exp_data.copy()

        if gas_phase:
            gas_gauss = self._get_gas_phase(x_range=x_range,
                                            exp_data=exp_data,
                                            gas_interval=gas_interval)
            df[self._GAS_PHASE] = gas_gauss
            df[self._CLEAN_EXP] -= df[self._GAS_PHASE]

        if contaminate or decontaminate:
            warnings.warn(UserWarning('Contaminants not implemented.'))

        if deconvolute:
            # Get a dict {specie: conc}, like sim_concs
            decon_concs = self._get_deconvolution(
                x_range=x_range,
                species=decon_species,
                species_guesses=sim_concs,
                exp_to_fit=df[self._CLEAN_EXP]
            )
            for dec_col in decon_cols:
                # dec_col is a tuple (self._EXPERIMENTAL, specie)
                df[dec_col] = self._get_gaussian(x_range=x_range,
                                                 specie=dec_col[1],
                                                 conc=decon_concs[dec_col[1]])
            df[self._DECON_ENV] = df[decon_cols].sum(axis=1)

        if simulate:
            for sim_col in sim_cols:
                # sim_col is a tuple (self._SIMULATED, specie)
                df[sim_col] = self._get_gaussian(x_range=x_range,
                                                 specie=sim_col[1],
                                                 conc=sim_concs[sim_col[1]])

            # Scales the gaussians and envelope by scale_factor
            # If self.autoscale, scale it to to self._get_autoscale()
            if autoscale:
                scale_factor = self._get_autoscale(df=df,
                                                   columns=sim_cols,
                                                   exp_data=exp_data)
                _echo.echo(f'Auto-scaling data to {scale_factor}...')
            for sim_col in sim_cols:
                df[sim_col] *= scale_factor

            df[self._ENVELOPE] = df[sim_cols].sum(axis=1)

        return df

    def _get_autoscale(self, df: pd.DataFrame,
                       columns: List[Tuple[str, sym.Symbol]],
                       exp_data: Union[pd.Series, None]) -> float:
        """Gets the factor by which to automatically scale.

        Currently, gives 1.0 unless there's an experimental, in which case it
        makes the peak experimental equal the peak envelope.

        Args:
            df: The pd.DataFrame with the gaussians to scale.
            columns: The columns in df you're going to scale.
            exp_data: Experimental data, if any, to autoscale to.
        """
        scale_factor = 1.0
        if exp_data is not None:
            envelope = df[columns].sum(axis=1)
            raw_max = max(envelope)
            scale_factor = max(exp_data) / raw_max

        return scale_factor

    def _get_deconvolution(self, x_range: np.ndarray,
                           species: List[sym.Symbol],
                           species_guesses: Dict[sym.Symbol, float],
                           exp_to_fit: pd.Series) -> Dict[sym.Symbol, float]:
        """Get the deconvolution of experimental into species.

        Args:
            x_range: np.ndarray, the x-range on which to calculate gaussians.
            species: A list of sym.Symbols, the species for which we're 
                guessing concentrations.
            species_guesses: A dict {species: conc} of initial guesses.
                This should be sim_concs, as far as I can tell.
            exp_to_fit: pd.Series, the experimental data to fit to.

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
                envelope += self._get_gaussian(x_range, specie, concs[index])
            return envelope

        fitted_concs, covariance = optimize.curve_fit(f=envelope_from_concs,
                                                      xdata=x_range,
                                                      ydata=exp_to_fit,
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
              autoscale: Optional[bool] = None,
              experimental: Optional[bool] = None,
              gas_phase: Optional[bool] = None,
              decontaminate: Optional[bool] = None,
              contaminate: Optional[bool] = None,
              deconvolute: Optional[bool] = None,
              # Args passed only to plotter
              sim_gaussians: bool = True,
              envelope: bool = True,
              experimental_raw: bool = True,
              experimental_clean: bool = True,
              decon_gaussians: bool = True,
              decon_envelope: bool = True,
              **kwargs) -> plt.Axes:
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
        if decontaminate is None:
            decontaminate = True
        if deconvolute is None:
            deconvolute = True

        # XPSExperiment.plot is overridden by Experiment.plot, but
        # this is an XPSObservable, so this doesn't make an infinite loop.
        return xps_obs.plot(ax=ax,
                            simulate=simulate,
                            sim_gaussians=sim_gaussians,
                            envelope=envelope,
                            experimental=experimental,
                            experimental_raw=experimental_raw,
                            experimental_clean=experimental_clean,
                            gas_phase=gas_phase,
                            deconvolute=deconvolute,
                            decon_gaussians=decon_gaussians,
                            decon_envelope=decon_envelope,
                            # TODO: Plot contaminants
                            **kwargs)

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
                 species: List[sym.Symbol] = None,
                 ignore: List[sym.Symbol] = None,
                 autoresample: bool = True,
                 autoscale: bool = True,
                 experimental: pd.Series = None,
                 gas_interval: Tuple[float, float] = None,
                 scale_factor: float = None,
                 title: str = '',
                 **options) -> XPSExperiment:
    """Simulate the given reaction system over time.

    Args:
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
        sol_t, sol_y = bulk_crn.solve_rsys_ode_till_eq(rsys, time, **options)
    else:
        sol_t, sol_y = bulk_crn.solve_rsys_ode(rsys, time, **options)

    cts = time_series.CRNTimeSeries(sol_t, sol_y, rsys)
    return cts.xps_with(species=species,
                        ignore=ignore,
                        autoresample=autoresample,
                        autoscale=autoscale,
                        experimental=experimental,
                        gas_interval=gas_interval,
                        scale_factor=scale_factor,
                        title=title)
