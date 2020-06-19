"""Utilities for manipulating XPS data.

Exports:
    XPSObservable: a container for simulated and experimental xps observables.

Example:
    Say you have a normal setup with a species manager, species x, y, and some
    experimental data in a pd.Series. The following is preferred:

    xps = XPSObservable({x: 2.1, y: 3.4}, sm, experimental=experimental_data,
                        gas_interval=(-1, 0), scale_factor=0.1)
    xps.plot()

    but is equivalent to the following, which resamples with each assignment:

    xps = XPSObservable({x: 2.1, y: 3.4}, sm)
    xps.experimental = experimental_data
    xps.gas_interval = (-1, 0)
    xps.scale_factor = 0.1
    xps.plot()
"""

from typing import Dict, List, Tuple, Union

import copy
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy import integrate
from scipy import stats
from sklearn import metrics
import sympy as sym

from lblcrn.crn_sym import species


class XPSObservable:
    """A container for a simulated observable of an XPS experiment.

    Attributes:
        df: The pandas DataFrame in which the observables are stored.
        species_concs: The concentrations of each species, for the creation of
            simulated data.
        species_manager: The SpeciesManager in use.
        title: A str used for the title during plotting.
    """

    _COLORS = ['red', 'green', 'orange', 'blue', 'purple', 'pink', 'yellow',
               'gray', 'cyan']
    _PLOT_MARGIN = 5
    _PLOT_RESOLUTION = 0.001
    _RESERVED_COLUMNS = ['envelope', 'experimental', 'gas_phase']
    _SIGMA = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)

    def __init__(self, species_concs: Dict[sym.Symbol, float],
                 species_manager: species.SpeciesManager,
                 experimental: pd.Series = None,
                 gas_interval: Tuple[float, float] = None,
                 scale_factor: float = 0,
                 title: str = ''):
        """Initialze the XPSObserable.

        Pass arguments instead of setting to prevent excessive re-sampling.

        Args:
            species_concs: A dictionary {specie: concentration} of the species
                and their concentrations.
            species_manager: A SpeciesManager, for binding energies.
            experimental: Optional, the experimental data.
            gas_interval: Optional, the interval in which to find the peak of
                the gas phase's gaussian.
            scale_factor: Optional, the factor by which to scale the gaussians.
            title: The title to name the default plot.
        """
        self.df = None
        self.species_concs = species_concs
        self.species_manager = species_manager
        self.title = title

        self._experimental = None
        if experimental is not None:
            self._experimental = experimental.copy()
        self._gas_interval = gas_interval
        self._scale_factor = scale_factor

        self.resample()

    # --- Accessors ----------------------------------------------------------

    @property
    def envelope(self):
        """The simulated envelope, the sum of all species' gaussians."""
        return self.df.envelope

    @property
    def experimental(self) -> Union[pd.Series, None]:
        """The experimental data. Might be None."""
        return self._experimental

    @experimental.setter
    def experimental(self, experimental: pd.Series):
        """Set the experimental data, prompting a resample."""
        self._experimental = experimental
        self.resample()

    @experimental.deleter
    def experimental(self):
        """Deletes experimental, prompting a resample."""
        self._experimental = None
        del self.df['experimental']
        del self.df['gas_phase']
        self.resample()

    @property
    def gas_interval(self) -> Union[Tuple[float, float], None]:
        """The interval in which the gas phase peak should be. May be None."""
        return self._gas_interval

    @gas_interval.setter
    def gas_interval(self, interval: Tuple[float, float]):
        """Sets the interval in which the gas phase peak should be.

        Give bounds for where the _peak_ is: if your bounds are too broad, it
        will get greedy, as it assumes that that interval is dominated by the
        gas phase. You can make lower equal to upper if you know the peak.

        Args:
            interval: a 2-tuple, with interval[0] <= interval[1].
        """
        if len(interval) != 2 or interval[0] > interval[1]:
            raise ValueError(f'Invalid interval {interval}')
        self._gas_interval = interval
        self.resample()

    @gas_interval.deleter
    def gas_interval(self):
        """Deletes the gas interval, prompting a resample."""
        self._gas_interval = None
        del self.df['gas_phase']
        self.resample()

    @property
    def gas_phase(self) -> Union[pd.Series, None]:
        """The gas phase part of the spectrum. Might be None."""
        if 'gas_phase' in self.df:
            return self.df.gas_phase
        else:
            return None

    @property
    def gaussians(self) -> pd.DataFrame:
        """The gaussians of the XPS observable.

        Doesn't include the envelope, experimental, or the gas phase.
        """
        gauss_cols = [col for col in self.df.columns
                      if col not in self._RESERVED_COLUMNS]
        return self.df[gauss_cols]

    @property
    def scale_factor(self) -> float:
        """The factor by which to scale the simulated data."""
        return self._scale_factor

    @scale_factor.setter
    def scale_factor(self, factor: float):
        """Sets the scale factor, prompting a resample."""
        self._scale_factor = factor
        self.resample()

    @scale_factor.deleter
    def scale_factor(self):
        """Deletes the scale factor, prompting a resample."""
        self._scale_factor = 0
        self.resample()

    @property
    def species(self) -> List[sym.Symbol]:
        """Give the species' symbols in this solution."""
        return list(self.species_concs.keys())

    @property
    def x_range(self) -> np.ndarray:
        """The x-values, energies, on which there is data."""
        return np.asarray(self.df.index)

    # --- Calculations -------------------------------------------------------

    def resample(self):
        """Recalculates the dataframe, in case anything updated.

        This is the only method which mutates self.df. It decides what to do
        based on:
        - if there's a experimental defined.
        - if there's a gas_interval defined.
        """
        x_range = self._get_x_range()
        self.df = pd.DataFrame(data=0, index=x_range, columns=['envelope'])

        # Add the experimental data and gas phase if possible.
        if self._experimental is not None:
            self.df['experimental'] = self._experimental

            gas_phase = self._get_gas_phase()
            if gas_phase is not None:
                self.df['gas_phase'] = gas_phase

        # Make the gaussians on the new x-range
        for specie in self.species:
            self.df[specie] = self._get_gaussian(specie)

        # Scale the gaussians
        if self._scale_factor <= 0:
            self._scale_factor = self._get_autoscale()
        for specie in self.species:
            self.df[specie] *= self._scale_factor

        # Add the scaled envelope
        self.df['envelope'] = self.gaussians.sum(axis=1)

    def _get_autoscale(self) -> float:
        """Gets the factor by which to automatically scale.

        Currently, gives one unless there's an experimental, in which case it
        makes the peak experimental equal the peak envelope."""
        if self.experimental is not None:
            raw_max = max(self.envelope) / self.scale_factor
            return max(self.experimental) / raw_max
        return 1.0

    def _get_gas_phase(self) -> Union[np.ndarray, None]:
        """Calculates the gas phase given the current experimental.

        Returns:
            None, if there isn't a gas interval or experimental, or an ndarray
            of the gas phase gaussian.
        """
        # If there is no gas interval or experimental, do nothing.
        if not self._gas_interval or self._experimental is None:
            return

        # Get the range to search for the gas peak in.
        search = self.experimental[self._gas_interval[0]:self._gas_interval[1]]

        # Get the location of the highest part of the experimental data in range
        peak_x = max(search.index, key=lambda index: search[index])
        peak = self.experimental[peak_x]

        # Make a gaussian the same height as the experimental gas phase peak
        gas_gaussian = stats.norm.pdf(self.x_range, peak_x, self._SIGMA)
        gas_gaussian *= (peak / max(gas_gaussian))

        return gas_gaussian

    def _get_gaussian(self, specie: sym.Symbol) -> np.ndarray:
        """Calculates the (unscaled) gaussian for the given species."""
        gaussian = np.zeros(self.x_range.size)

        for orbital in self.species_manager[specie].orbitals:
            gaussian += self.species_concs[specie] * orbital.splitting * \
                stats.norm.pdf(self.x_range, orbital.binding_energy,
                               self._SIGMA)

        return gaussian

    def _get_x_range(self) -> np.ndarray:
        """Gets an adequate x-range on which to calculate gaussians.

        If there's experimental data, it'll just return the x-range of it. If
        there isn't, then it will pick a range which contains all the binding
        energies in question.

        Returns:
            An ndarray on which we're going to calculate gaussians.
        """
        # If there's experimental data, we have to use that x-range.
        if self._experimental is not None:
            return np.asarray(self._experimental.index)

        # Otherwise, pick intelligently based on binding energies.
        binding_energies = []
        for specie in self.species:
            for orbital in self.species_manager[specie].orbitals:
                binding_energies.append(orbital.binding_energy)

        x_lower = min(binding_energies) - self._PLOT_MARGIN
        x_upper = max(binding_energies) + self._PLOT_MARGIN
        x_range = np.arange(x_lower, x_upper, self._PLOT_RESOLUTION)

        return x_range

    # --- Plotting -----------------------------------------------------------

    def plot(self, show=True):
        """Plot the XPS observable."""
        for index, specie in enumerate(self.gaussians):
            plt.fill(self.x_range, self.gaussians[specie], label=specie,
                     color=self._COLORS[index])

        if self.gas_phase is not None:
            plt.fill(self.x_range, self.gas_phase, label='gas phase',
                     color='gray')

        plt.plot(self.x_range, self.envelope, color='black', linewidth=4)

        if self.experimental is not None:
            plt.plot(self.x_range, self.experimental, color='green')

        plt.legend()
        plt.title(self.title)
        plt.gca().invert_xaxis()  # XPS Plots are backwards

        if show:
            plt.show()

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

    # --- Utility ------------------------------------------------------------

    def _repr_html(self):
        """For iPython, mostly just gives the DataFrame."""
        # TODO: do this without accessing private member?
        # TODO: Also it doesn't appera to get called.
        return self.df._repr_html_()