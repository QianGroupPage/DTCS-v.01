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
               'gray', 'cyan']
    _RESERVED_COLUMNS = ['envelope', 'experimental', 'gas_phase']

    def __init__(self, df: pd.DataFrame, title: str = ''):
        self.df = df
        self.title = title

    @property
    def envelope(self) -> pd.Series:
        """The simulated envelope, the sum of all species' gaussians."""
        return self.df.envelope

    @property
    def experimental(self) -> Optional[pd.Series]:
        """The experimental data. Might be None."""
        if 'experimental' in self.df:
            return self.df.experimental
        else:
            return None

    @property
    def gas_phase(self) -> Optional[pd.Series]:
        """The gas phase part of the spectrum. Might be None."""
        if 'gas_phase' in self.df:
            return self.df.gas_phase
        else:
            return None

    @property
    def gaussians(self) -> pd.DataFrame:
        """The gaussians of the XPS observable."""
        gauss_cols = [col for col in self.df.columns
                      if col not in self._RESERVED_COLUMNS]
        return self.df[gauss_cols]

    @property
    def x_range(self) -> np.ndarray:
        """The x-values, energies, on which there is data."""
        return np.asarray(self.df.index)

    def plot(self, ax: plt.Axes, **kwargs):
        """Plot the XPS observable.

        Args:
            ax: The plt.Axes on which to plot.
            **kwargs: Forwarded.
        """
        # Sort the Gaussians before plotting to overlay smaller peaks on top of larger ones
        for index, specie in sorted(enumerate(self.gaussians), key=lambda x: max(self.gaussians[x[1]]), reverse=True):
            ax.fill(self.x_range, self.gaussians[specie], label=specie,
                    color=self._COLORS[index])

        if self.gas_phase is not None:
            ax.fill(self.x_range, self.gas_phase, label='gas phase',
                    color='gray')

        ax.plot(self.x_range, self.envelope, color='black', linewidth=4)

        if self.experimental is not None:
            ax.plot(self.x_range, self.experimental, color='green')

        ax.legend()
        ax.set_title(self.title)
        ax.invert_xaxis()  # XPS Plots are backwards

    # --- Data Analysis ------------------------------------------------------

    def rmse(self):
        return np.sqrt(metrics.mean_squared_error(self.experimental,
                                                  self.envelope))

    def mae(self):
        return metrics.mean_absolute_error(self.experimental, self.envelope)

    def integral_diff_outside_experimental(self):
        return integrate.trapz(self.envelope - self.experimental, self.x_range)

    def integral_diff_inside_experimental(self):
        return integrate.trapz(self.experimental - self.envelope, self.x_range)

    def integral_diff_between(self):
        return integrate.trapz(np.absolute(self.experimental - self.envelope),
                               self.x_range)

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
                 species_concs: Dict[sym.Symbol, float],
                 species_manager: species.SpeciesManager,
                 autoresample: bool = True,
                 autoscale: bool = True,
                 experimental: pd.Series = None,
                 gas_interval: Tuple[float, float] = None,
                 scale_factor: float = 1.0,
                 title: str = ''):
        """Initialze the XPSExperiment.

        Pass arguments instead of setting to prevent excessive re-sampling.

        Args:
            species_concs: A dictionary {specie: concentration} of the species
                and their concentrations.
            species_manager: A SpeciesManager, for binding energies.
            autoresample: Defaults to true, decides if it resamples on edits.
            autoscale: Defaults to true, decides if it will automatically scale
                the gaussians and envelope to match the experimental data.
            experimental: Optional, the experimental data.
            gas_interval: Optional, the interval in which to find the peak of
                the gas phase's gaussian.
            scale_factor: Optional, the factor by which to scale the gaussians.
            title: The title to name the default plot.
        """
        # This doesn't call XPSObservable.__init__(), which is correct.
        super().__init__()
        self.species_concs = species_concs
        self.species_manager = species_manager
        self.autoresample = autoresample
        self.autoscale = autoscale
        self.title = title

        # Experimental data-related
        self._experimental = None
        if experimental is not None:
            self._experimental = experimental.copy()
        self._gas_interval = gas_interval

        # Scaling-related
        self._scale_factor = scale_factor
        if scale_factor != 0.0:
            self.autoscale = False

        self._autoresample()

    # --- Accessors ----------------------------------------------------------

    @XPSObservable.experimental.setter
    def experimental(self, experimental: pd.Series):
        """Forwards to self.set_experimental."""
        self.set_experimental(experimental)

    @experimental.deleter
    def experimental(self):
        """Forwards to self.del_experimental."""
        self.del_experimental()

    @property
    def gas_interval(self) -> Optional[Tuple[float, float]]:
        """The interval in which the gas phase peak should be. May be None."""
        return self._gas_interval

    @gas_interval.setter
    def gas_interval(self, interval: Tuple[float, float]):
        """Forwards to self.set_gas_interval."""
        if len(interval) != 2:
            raise ValueError(f'Invalid interval {interval}')
        self.set_gas_interval(lower=interval[0], upper=interval[1])

    @gas_interval.deleter
    def gas_interval(self):
        """Forwards to self.del_gas_interval"""
        self.del_gas_interval()

    @property
    def scale_factor(self) -> float:
        """The factor by which to scale the simulated data."""
        return self._scale_factor

    @scale_factor.setter
    def scale_factor(self, factor: float):
        """Forwards to self.set_scale_factor."""
        self.set_scale_factor(factor)

    @scale_factor.deleter
    def scale_factor(self):
        """Forwards to self.del_scale_factor."""
        self.del_scale_factor()

    @property
    def species(self) -> List[sym.Symbol]:
        """Give the species' symbols in this solution."""
        return list(self.species_concs.keys())

    # --- Mutators -----------------------------------------------------------
    # If autoresample is on, these will prompt an overwriting resample.

    def set_experimental(self, experimental: pd.Series):
        """Set the experimental data, prompting an autoresample."""
        self._experimental = experimental
        self._autoresample()

    def del_experimental(self):
        """Deletes experimental, prompting an autoresample.."""
        self._experimental = None
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

    def _autoresample(self):
        """Does an overwriting resample, if necessary."""
        if self.autoresample:
            _echo.echo('Auto-resampling data...')
            self.resample(overwrite=True)

    # --- Calculations -------------------------------------------------------
    # Note that the only function in this section which modifies the state is
    # self.resample(), and that's only if overwrite=True.

    def resample(self, overwrite=True, species=None,
                 ignore=[]) -> XPSObservable:
        """Recalculates the dataframe in case anything updated.

        Args:
            overwrite: If true, overwrites self.df with the resampled df.
            species: A list of sym.Symbols, if you only want to sample for
                those species.
            ignore: A list of sym.Symbols to not include.

        Returns:
            An XPSObservable with the resampled data.
        """  # TODO(Andrew) Typehints?
        species = self._get_species_not_ignored(species, ignore + self.ignore)

        x_range = self._get_x_range(species)
        df = pd.DataFrame(data=0, index=x_range, columns=['envelope'])

        # Add the experimental data and gas phase if possible.
        if self._experimental is not None:
            df['experimental'] = self._experimental

            gas_phase = self._get_gas_phase(x_range)
            if gas_phase is not None:
                df['gas_phase'] = gas_phase

        # Make the gaussians on the new x-range
        for specie in species:
            df[specie] = self._get_gaussian(specie, x_range)

        # Scales the gaussians and envelope by self._scale_factor
        # If self.autoscale, scale it to to self._get_autoscale()
        scale = self._scale_factor
        if self.autoscale:
            scale = self._get_autoscale(df, species)
        for specie in species:
            df[specie] *= scale

        df['envelope'] = df[species].sum(axis=1)

        if overwrite:
            self.df = df
        return XPSObservable(df)

    def _get_autoscale(self, df: pd.DataFrame,
                       species: List[sym.Symbol]) -> float:
        """Gets the factor by which to automatically scale.

        Currently, gives 1.0 unless there's an experimental, in which case it
        makes the peak experimental equal the peak envelope.

        Args:
            df: The pd.DataFrame with the gaussians to scale.
            species: The columns of df which are scale-able gaussians.
        """
        scale = 1.0
        if self._experimental is not None:
            envelope = df[species].sum(axis=1)
            raw_max = max(envelope)
            scale = max(self._experimental) / raw_max

        _echo.echo(f'Auto-scaling data to {scale}...')
        return scale

    def _get_gas_phase(self, x_range: np.ndarray) -> Optional[np.ndarray]:
        """Calculates the gas phase given the current experimental.

        Args:
            x_range: The x-values on which to caluclate the gaussian.

        Returns:
            None, if there isn't a gas interval or experimental, or an ndarray
            of the gas phase gaussian.
        """
        # If there is no gas interval or experimental, do nothing.
        if not self._gas_interval or self._experimental is None:
            return

        # Get the range to search for the gas peak in.
        search = self._experimental[self._gas_interval[0]:self._gas_interval[1]]

        # Get the location of the highest part of the experimental data in range
        peak_x = max(search.index, key=lambda index: search[index])
        peak = self._experimental[peak_x]

        # Make a gaussian the same height as the experimental gas phase peak
        gas_gaussian = stats.norm.pdf(x_range, peak_x, self._SIGMA)
        gas_gaussian *= (peak / max(gas_gaussian))

        return gas_gaussian

    def _get_gaussian(self, specie: sym.Symbol,
                      x_range: np.ndarray) -> np.ndarray:
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
            gaussian += self.species_concs[specie] * orbital.splitting * \
                        stats.norm.pdf(x_range, orbital.binding_energy,
                                       self._SIGMA)

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
        if self._experimental is not None:
            return np.asarray(self._experimental.index)

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

    def _plot(self, species: List[sym.Symbol], ax: plt.Axes, **kwargs):
        """Plot the XPS observable.

        Args:
            ax: The plt.Axes on which to plot.
            species: A list of sym.Symbols, the species to plot.
            **kwargs: Forwarded.
        """
        xps_obs = self.resample(overwrite=False, species=species)
        xps_obs.plot(ax=ax, **kwargs)

    # --- Utility -------------------------------------------------------------

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict()
        d['species_concs'] = {str(symbol): conc for symbol, conc in
                              self.species_concs.items()}
        d['species_manager'] = self.species_manager.as_dict()
        d['autoresample'] = self.autoresample
        d['autoscale'] = self.autoscale
        d['experimental'] = self._experimental.to_json() if \
            self._experimental is not None else None
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


def simulate_xps(rsys: reaction.RxnSystem, time: float = 1,
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
    # TODO(Andrew): Solve at equilibrium when no time is specified.
    sol_t, sol_y = bulk_crn.solve_rsys_ode(rsys, time, **options)
    sol = time_series.CRNTimeSeries(sol_t, sol_y, rsys)
    return sol.xps_with(species=species,
                        ignore=ignore,
                        autoresample=autoresample,
                        autoscale=autoscale,
                        experimental=experimental,
                        gas_interval=gas_interval,
                        scale_factor=scale_factor,
                        title=title)
