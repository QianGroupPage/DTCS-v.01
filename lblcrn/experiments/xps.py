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
    _DECONV_ENV = (_EXPERIMENTAL, 'deconv_envelope')
    _CLEAN_EXP = (_EXPERIMENTAL, 'clean')
    _RAW = (_EXPERIMENTAL, 'raw')

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
        if self._RAW in self.df:
            return self.df[self._RAW]
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
        if self._DECONV_ENV in self.df:
            return self.df[self._DECONV_ENV]
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
             simulated=True,
             experimental=True,
             sim_gaussians=True,
             deconvoluted=True,
             deconv_gaussians=True,
             deconv_envelope=True,
             gas_phase=True,
             envelope=True,
             experimental_clean=True,
             experimental_raw=True, # TODO: perhaps rename these all to show_ or show=?
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
        if not simulated:
            sim_gaussians = False
            envelope = False
        if not experimental:
            deconvoluted = False
            gas_phase = False
            experimental_raw = False
            experimental_clean = False
        if not deconvoluted:
            deconv_envelope = False
            deconv_gaussians = False

        if sim_gaussians and self.gaussians_simulated is not None:
            for index, column in enumerate(sorted(self.gaussians_simulated,
                                                  key=gauss_col_sort_key,
                                                  reverse=True)):
                ax.fill(self.x_range, self.df[column],
                        label=f'Sim. {column[1]}',
                        color=self._COLORS[index])

        if deconv_gaussians and self.gaussians_deconvoluted is not None:
            for index, column in enumerate(sorted(self.gaussians_deconvoluted,
                                                  key=gauss_col_sort_key,
                                                  reverse=True)):
                ax.plot(self.x_range, self.df[column],
                        label=f'Deconv. {column[1]}',
                        color=self._COLORS[index])

        if deconv_envelope and self.decon_envelope is not None:
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
        species_concs: The concentrations of each species, for the creation of
            simulated data.
        species_manager: The SpeciesManager in use.
    """

    _PLOT_MARGIN = 5
    _PLOT_RESOLUTION = 0.001
    _SIGMA = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)

    def __init__(self,
                 species_manager: species.SpeciesManager,
                 title: str = '',
                 species: List[sym.Symbol] = None,
                 species_concs: Dict[sym.Symbol, float] = None,
                 scale_factor: float = 0.0,
                 experimental: pd.Series = None,
                 gas_interval: Tuple[float, float] = None,
                 autoresample: bool = True,
                 autoscale: bool = True,):  # TODO: Add init resample options
        """Initialze the XPSExperiment.

        Pass arguments instead of setting to prevent excessive re-sampling.

        Args:
            species_manager: A SpeciesManager, for binding energies.
            species_concs: A dictionary {specie: concentration} of the species
                and their concentrations. This or experimental is required.
            experimental: The experimental data. This or species_concs is
                required.
            autoresample: Defaults to true, decides if it resamples on edits.
            autoscale: Defaults to true, decides if it will automatically scale
                the gaussians and envelope to match the experimental data.
            gas_interval: Optional, the interval in which to find the peak of
                the gas phase's gaussian.
            scale_factor: Optional, the factor by which to scale the gaussians.
            title: The title to name the default plot.
        """
        super().__init__()
        XPSObservable.__init__(self, species_manager=species_manager,
                               df=None, title=title)

        # Deal with required arguments
        if not (species_concs or experimental is not None):
            raise ValueError(f'{self.__class__.__name__} needs at least'
                             f'species_concs or experimental defined.')

        if not species_concs:
            species_concs = {}
        if experimental is not None:
            experimental = experimental.copy()
        if not species:
            if species_concs:
                species = list(species_concs.keys())
            else:
                # TODO: perhaps an echo saying that we're assuming this?
                species = species_manager.species

        # Default scale factor to 1
        # If they pick 1.0 on their own, then it disables autoscale.
        if scale_factor == 0.0:
            scale_factor = 1.0
        else:
            autoscale = False

        # Flags
        self.autoresample = autoresample
        self.autoscale = autoscale
        # Internal data, meant to act like a descriptor of state.
        self._species = species
        self.species_concs = species_concs
        self._exp_data = experimental
        self._gas_interval = gas_interval
        self._scale_factor = scale_factor

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

    def set_scale_factor(self, factor: float):
        """Sets the scale factor, disabling autoscale and prompting an
        autoresample."""
        self.autoscale = False
        self._scale_factor = factor
        self._autoresample()

    def del_scale_factor(self):
        """Deletes the scale factor, enabling autoscale and prompting an
        autoresample."""
        self.autoscale = True
        self._scale_factor = 0
        self._autoresample()

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
            self._species.pop(specie)
        self._autoresample()

    def _autoresample(self):
        """Does an overwriting resample, if necessary."""
        if self.autoresample:
            _echo.echo('Auto-resampling data...')
            self.resample(overwrite=True)

    # --- Calculations -------------------------------------------------------
    # Note that functions in this section will only modify state if
    # overwwrite=True.

    def resample(self, overwrite=True, species=None, ignore=None,
                 experimental=True,
                 simulated=True,
                 gas_phase=True,
                 deconvolute=True) -> XPSObservable:  # TODO(Andrew) Typehints?
        """Recalculates the dataframe in case anything updated.

        Args:
            overwrite: If true, overwrites self.df with the resampled df.
            species: A list of sym.Symbols, if you only want to sample for
                those species.
            ignore: A list of sym.Symbols to not include.
            TODO: add rest of params

        Returns:
            An XPSObservable with the resampled data.
        """
        species = self._get_species_not_ignored(species, ignore)

        # This uses the experimental data, even if experimental=False
        # I'm leaving it in because I can't think of a time you would want
        # to have the x-axes mismatch.
        x_range = self._get_x_range(species)
        scale = self._scale_factor

        # TODO: Warn the user if I override what they asked
        if not self.species_concs:
            simulated = False
        if self._exp_data is None:
            experimental = False
        if not self._gas_interval or not experimental:
            gas_phase = False
        if not experimental or not species:
            deconvolute = False

        # Make the pd.Index-s
        columns = []
        if experimental:
            columns.append(self._RAW)
            if gas_phase:
                columns.append(self._GAS_PHASE)
            if deconvolute:
                dec_cols = [(self._EXPERIMENTAL, specie) for specie in species]
                columns.extend(dec_cols)
        if simulated:
            columns.append(self._ENVELOPE)
            sim_cols = [(self._SIMULATED, specie) for specie in species]
            columns.extend(sim_cols)

        row_index = pd.Index(x_range, name='eV')
        col_index = pd.MultiIndex.from_tuples(columns)

        df = pd.DataFrame(data=0, index=row_index, columns=col_index)

        # Add the experimental data and gas phase if possible.
        if experimental:
            df[self._RAW] = self._exp_data.copy()
            df[self._CLEAN_EXP] = self._exp_data.copy()

        if gas_phase:
            gas_gauss = self._get_gas_phase(x_range)
            df[self._GAS_PHASE] = gas_gauss
            df[self._CLEAN_EXP] -= df[self._GAS_PHASE]

        if deconvolute:
            noise = np.zeros(x_range.size)
            if gas_phase:
                noise += df[self._GAS_PHASE]

            # A dict {specie: conc}, like self.species_concs
            decon_concs = self._get_deconvolution(x_range, species, noise)
            for dec_col in dec_cols:
                # dec_col is a tuple (self._EXPERIMENTAL, specie)
                df[dec_col] = self._get_gaussian(x_range,
                                                 dec_col[1],
                                                 decon_concs[dec_col[1]])

            df[self._DECONV_ENV] = df[dec_cols].sum(axis=1)

        # Make the gaussians on the new x-range
        if simulated:
            for sim_col in sim_cols:
                # sim_col is a tuple (self._SIMULATED, specie)
                df[sim_col] = self._get_gaussian(x_range,
                                                 sim_col[1],
                                                 self.species_concs[sim_col[1]])

            # Scales the gaussians and envelope by self._scale_factor
            # If self.autoscale, scale it to to self._get_autoscale()
            if self.autoscale:
                scale = self._get_autoscale(df, sim_cols, experimental)
            for sim_col in sim_cols:
                df[sim_col] *= scale

            df[self._ENVELOPE] = df[sim_cols].sum(axis=1)

        if overwrite:
            self._scale_factor = scale
            self.df = df
        return XPSObservable(species_manager=self.species_manager,
                             df=df,
                             title=self.title)

    def _get_autoscale(self, df: pd.DataFrame, columns: Tuple[str, sym.Symbol],
                       experimental=False) -> float:
        """Gets the factor by which to automatically scale.

        Currently, gives 1.0 unless there's an experimental, in which case it
        makes the peak experimental equal the peak envelope.

        Args:
            df: The pd.DataFrame with the gaussians to scale.
            species: The columns of df which are scale-able gaussians.
        """
        scale = 1.0
        if experimental:
            envelope = df[columns].sum(axis=1)
            raw_max = max(envelope)
            scale = max(self._exp_data) / raw_max

        _echo.echo(f'Auto-scaling data to {scale}...')
        return scale

    def _get_deconvolution(self, x_range: np.ndarray,
                           species: List[sym.Symbol],
                           noise: pd.Series) -> Dict[sym.Symbol, float]:
        """Get the deconvolution of experimental into species.

        Args:
            species: The species for which we're guessing concentrations.
            noise: Stuff to exclude from the experimental data.
                e.g. the gas phase, contaminants, etc.

        Returns:
            A dictionary of {species: conc} of the guessed concentration.
        """
        to_fit = self._exp_data - noise

        # Guess 1 for everything, unless there's simulated data, then use that.
        conc_guesses = np.ones(len(species))
        for index, specie in enumerate(species):
            if specie in self.species_concs:
                conc_guesses[index] = self.species_concs[specie]
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
                                                      ydata=to_fit,
                                                      p0=conc_guesses,
                                                      bounds=(0, np.inf))
        return dict(zip(species, fitted_concs))


    def _get_gas_phase(self, x_range: np.ndarray) -> np.ndarray:
        """Calculates the gas phase given the current experimental.

        Args:
            x_range: The x-values on which to caluclate the gaussian.

        Returns:
            An ndarray of the gas phase gaussian.
        """
        # Get the range to search for the gas peak in.
        search = self._exp_data[self._gas_interval[0]:self._gas_interval[1]]

        # Get the location of the highest part of the experimental data in range
        peak_x = max(search.index, key=lambda index: search[index])
        peak = self._exp_data[peak_x]

        # Make a gaussian the same height as the experimental gas phase peak
        gas_gaussian = stats.norm.pdf(x_range, peak_x, self._SIGMA)
        gas_gaussian *= (peak / max(gas_gaussian))

        return gas_gaussian

    def _get_gaussian(self, x_range: np.ndarray, specie: sym.Symbol,
                      conc: float) -> np.ndarray:
        """Calculates the (unscaled) gaussian for the given species.

        Args:
            specie: The specie for which to get the gaussian.
            x_range: The x-values on which to caluclate the gaussian.

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
        """Gets an adequate x-range on which to calculate gaussians.

        If there's experimental data, it'll just return the x-range of it. If
        there isn't, then it will pick a range which contains all the binding
        energies in question.

        Args:
            species: The species to include in the x_range.

        Returns:
            An ndarray on which we're going to calculate gaussians.
        """
        # If there's experimental data, we have to use that x-range.
        if self._exp_data is not None:
            return np.asarray(self._exp_data.index)

        # Otherwise, pick intelligently based on binding energies.
        binding_energies = []
        for specie in species:
            for orbital in self.species_manager[specie].orbitals:
                binding_energies.append(orbital.binding_energy)

        x_lower = min(binding_energies) - self._PLOT_MARGIN
        x_upper = max(binding_energies) + self._PLOT_MARGIN
        x_range = np.arange(x_lower, x_upper, self._PLOT_RESOLUTION)

        return x_range

    # --- Plotting -----------------------------------------------------------

    def _plot(self, species: List[sym.Symbol], ax: plt.Axes,
              experimental=True, simulated=True,
              gas_phase=True, deconvolute=True, **kwargs) -> plt.Axes:
        """Plot the XPS observable.

        Args:
            ax: The plt.Axes on which to plot.
            species: A list of sym.Symbols, the species to plot.
            **kwargs: Forwarded.

        Returns:
            The plt.Axes of the plot.
        """
        # This makes an XPSObservable.
        xps_obs = self.resample(overwrite=False, species=species,
                                experimental=experimental,
                                simulated=simulated,
                                gas_phase=gas_phase,
                                deconvolute=deconvolute)
        # XPSExperiment.plot is overridden by Experiment.plot, but
        # this is an XPSObservable, so this doesn't make an infinite loop.
        return xps_obs.plot(ax=ax,
                            experimental=experimental,
                            simulated=simulated,
                            gas_phase=gas_phase,
                            deconvolute=deconvolute,
                            **kwargs)

    # --- Utility -------------------------------------------------------------

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict()
        d['species_concs'] = {str(symbol): conc for symbol, conc in
                              self.species_concs.items()}
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
                 scale_factor: float = 0.0,
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
