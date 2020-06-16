import bisect
import collections
from typing import Dict, List, Tuple, Union
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import integrate
from scipy import stats
from sklearn import metrics
import sympy as sym

_SIGMA = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)
_COLORS = ['red', 'green', 'orange', 'blue', 'purple', 'pink', 'yellow', 'gray', 'cyan']
_PLOT_MARGIN = 5
_PLOT_RESOLUTION = 0.001


class XPSObservable:
    """A container for a simulated observable of an XPS experiment.

    Attributes:
        title: A str used for the title during plotting.
    """

    _RESERVED_COLUMNS = ['envelope', 'experimental', 'gas_phase']

    def __init__(self, sim=None):
        """Initializes everything to none; won't work until you call update."""
        self.title = ''
        self._df = None
        self._gas_interval = ()
        self._sim = sim

    @property
    def df(self) -> pd.DataFrame:
        """The pandas DataFrame describing the observable.

        Columns 'envelope', 'experimental', and 'gas_phase' are not gaussians,
        but everything else is.

        Note: Accessing it might prompt the DataFrame to be re-evaluated.
        """
        self._check_resim()
        return self._df

    @property
    def envelope(self):
        """The simulated envelope, the sum of all species' gaussians.

        Note: Accessing it might prompt the DataFrame to be re-evaluated.
        """
        self._check_resim()
        return self._df.envelope

    @property
    def experimental(self) -> Union[pd.Series, None]:
        """The experimental data. Might be None."""
        if self._df is not None and 'experimental' in self._df:
            return self._df.experimental
        return None

    @experimental.setter
    def experimental(self, experimental: pd.Series):
        """Sets the experimental data, updating the simulated data if necessary.

        If you give a x-range or DataFrame with index different from that of
        the simulated data, it will drop the simulated data and throw a
        warning, resetting to zeroes on that x-range.

        Args:
            experimental: An ndarray of the experimental spectrum, with index.
        """
        x_range = np.asarray(experimental.index)

        # Reset self._df if it has never been initialized.
        if self._df is None:
            self._reset_df(x_range)

        # If the x-axes are irreconcilable, reset and throw a warning.
        if self.x_range.size != x_range.size or (self.x_range != x_range).any():
            self._reset_df(x_range)
            warnings.warn('x-ranges differ, dropped existing data.')

        self._df['experimental'] = experimental

        self._scale_gaussians()
        self._calc_gas_phase()

    @experimental.deleter
    def experimental(self):
        """Deletes experimental but also gas_phase."""
        del self._df['experimental']
        del self._df['gas_phase']

    @property
    def gas_phase(self) -> Union[pd.Series, None]:
        """The gas phase part of the spectrum. Might be None."""
        if self._df is not None and 'gas_phase' in self._df:
            return self._df.gas_phase
        return None

    @property
    def gas_interval(self) -> Tuple:
        """The interval in which the gas phase peak should be."""
        return self._gas_interval

    @gas_interval.setter
    def gas_interval(self, interval):
        """Sets the interval in which the gas phase peak should be.

        Give bounds for where the _peak_ is: if your bounds are too broad, it
        will get greedy, as it assumes that that interval is dominated by the
        gas phase. You can make lower equal to upper if you know the peak.

        Args:
            interval: a 2-tuple, with interval[0] <= interval[1].
        """
        if len(interval) != 2 or interval[0] > integrate[1]:
            raise ValueError(f'Invalid interval {interval}')
        self._gas_interval = interval
        self._calc_gas_phase()

    @property
    def gaussians(self) -> pd.DataFrame:
        """The gaussians of the XPS observable.

        Doesn't include the envelope, experimental, or the gas phase.

        Note: Accessing it might prompt the DataFrame to be re-evaluated.
        """
        self._check_resim()
        gauss_cols = [col for col in self._df.columns
                      if col not in self._RESERVED_COLUMNS]
        return self._df[gauss_cols]

    @gaussians.setter
    def gaussians(self, gaussians: pd.DataFrame):
        """Sets the simulated gaussians.

        This will prompt a calculation of the envelope. Don't use the names
        'envelope', 'experimental', or 'gas_phase' in your DataFrame.

        Args:
            gaussians: a DataFrame of all the gaussians.

        Raises:
            ValueError: if you give a DataFrame with an index which doesn't
                match the experimental data.
        """
        experimental = self.experimental
        self._df = gaussians
        if experimental:
            self._df['experimental'] = experimental

        self._scale_gaussians()
        self._calc_gas_phase()

    @property
    def x_range(self) -> np.ndarray:
        """The x-values, energies, on which there is data.

        Note: Accessing it might prompt the DataFrame to be re-evaluated.
        """
        if self._df is None:
            self._proper_resim()
        return np.asarray(self._df.index)

    def plot(self, species: List[sym.Symbol] = None,
             ignore: List[sym.Symbol] = None, t=-1):
        """ TODO

        Args:
            species:
            ignore:
            t:
        """
        # Handles parameters that sim might have to take.
        self._prompt_resim(species=species, ignore=ignore, t=t)

        # Plot everything which is defined.
        for index, specie in enumerate(self.gaussians):
            plt.fill(self.x_range, self.gaussians[specie], label=specie,
                     color=_COLORS[index])

        if self.gas_phase is not None:
            plt.fill(self.x_range, self.gas_phase, label='gas phase', color='gray')

        plt.plot(self.x_range, self.envelope, color='black', linewidth=4)

        if self.experimental is not None:
            plt.plot(self.x_range, self.experimental, color='green')

        plt.legend()
        plt.title(self.title)
        plt.gca().invert_xaxis()  # XPS Plots are backwards
        plt.show()

    def _calc_envelope(self):
        """Sets self._df['envelope'] to be the sum of the gaussians."""
        self._df['envelope'] = self.gaussians.sum(axis=1)

    def _calc_gas_phase(self):
        """Calculates the gas phase given the current experimental."""
        # If there is no gas interval or experimental, do nothing.
        if not (self._gas_interval and self.experimental is not None):
            return

        # Get the part of the x_range that the gas_interval corresponds to
        # i.e., x_range[lower_index] will be very close to gas_interval[0]
        # bisect is a standard library binary search function.
        lower_index = bisect.bisect_left(self.x_range, self._gas_interval[0])
        upper_index = bisect.bisect_left(self.x_range, self._gas_interval[1]) + 1

        # Get the location of the highest part of the experimental data in range
        peak_index = max(range(lower_index, upper_index),
                         key=lambda index: self.experimental[index])
        peak = self.experimental[peak_index]

        # Make a gaussian the same height as the experimental gas phase peak
        gas_gaussian = stats.norm.pdf(self.x_range, self.x_range[peak_index],
                                      _SIGMA)
        gas_gaussian *= (peak / max(gas_gaussian))

        self._df['gas_phase'] = gas_gaussian

    def _scale_gaussians(self):
        """Scale the gaussians so that the envelope and experimental have the
        same maximum. Also calculates envelope."""
        self._calc_envelope()

        # Prevent a division by zero error.
        if (self.envelope == 0).all():
            return

        # Scale it to the max of experimental, or 1.
        scale = 1 / max(self.envelope)
        if self.experimental:
            scale = max(self.experimental) / max(self.envelope)

        gauss_cols = [col for col in self._df.columns
                      if col not in self._RESERVED_COLUMNS]
        self._df[gauss_cols] = self._df[gauss_cols] * scale

        # Recalculate envelope to deal with new gaussians
        self._calc_envelope()

    def _check_resim(self, **kwargs):  # TODO
        if self._df is None or (self.envelope == 0).all():
            self._prompt_resim(**kwargs)

    def _prompt_resim(self, **kwargs):  # TODO
        if self._sim:
            self._sim.calc_xps(**kwargs)

    def _reset_df(self, x_range: np.ndarray):
        """Resets the dataframe to zeros to maintain invariants.

        Args:
            x_range: A float index, representing energies on the x-axis.
        """

        self._df = pd.DataFrame(data=0, index=x_range, columns=['envelope'])

class Solution:
    """A time series solution of an chemical reaction ODE, with utilities.

    Defines a solution to a system of differential equations, providing a
    variety of functions to manipulate and visualize the results

    Attributes:
        df: The pandas DataFrame associated with the time series.
        xps: An XPSObservable object. It will prompt the Solution to update it
            if you ask for something that doesn't exist yet.
    """
    
    def __init__(self, t: List[float], y: List[List[float]], rsys):

        self.df = pd.DataFrame(data=np.transpose(y),
                               index=pd.Index(t, name='time'),
                               columns=pd.Index(rsys.get_symbols(),
                                                name='species'))
        self.xps = XPSObservable(sim=self)

        self.rsys = rsys
        self.species_manager = rsys.species_manager

        # TODO ---------------------------------------------------------------

        self._default_ignore = []

        # The time range of the reaction
        self.t = t
        # A map of symbols to concentrations over time
        self.states: Dict[str, List[float]] = {}
        # A map of symbols to substance information
        self.substances = dict(zip(rsys.get_symbols(), rsys.get_species()))
        #self.substances: Dict[str, species.Species] = {}
        for symbol, sub in self.substances.items():
        #    self.substances[substances[i].symbol] = substances[i]
            self.states[symbol] = y[rsys.symbol_index[symbol]]

        self.envelope = []
        self.binding_energies = []
        self.resampled_intensity = []
        self.distributions = []
        #self.xps = None

    # --- Accessor Methods ---------------------------------------------------

    def at(self, t: float):
        """Gives the state of the system at time <= t, or at end if t < 0."""
        if t < 0:
            t = self.t[-1]
        return self.df.iloc[self._time_to_index(t)]

    @property
    def species(self) -> List[sym.Symbol]:
        """Give the species' symbols in this solution."""
        return self.rsys.get_symbols()

    # --- Calculating Observables --------------------------------------------

    def calc_xps(self, species: List[sym.Symbol] = None,
                 ignore: List[sym.Symbol] = None, t: float = -1):
        """Calculates a simulated XPS observable at time t.

        Doesn't return; instead, saves the information into self.xps.

        Args:
            species: The speices to include in the caluclation.
            t: The time at which to take a snapshot.
        """
        species = self._get_species_not_ignored(species, ignore)

        # Pick x_range: if there is experimental data, use that. Otherwise,
        # pick it so that it contains just the binding energies of the species
        # in question.
        if self.xps.experimental is not None:
            x_range = self.xps.x_range
        else:
            binding_energies = []
            for specie in species:
                binding_energies.extend(self._get_binding_energies(specie))

            x_lower = min(binding_energies) - _PLOT_MARGIN
            x_upper = max(binding_energies) + _PLOT_MARGIN
            x_range = np.arange(x_lower, x_upper, _PLOT_RESOLUTION)

        # Make a gaussian for each specie in question.
        gaussians = [self._get_species_gauss(specie, x_range, t) for specie
                     in self.species]
        names = [self.species_manager[specie].name for specie in self.species]

        xps_df = pd.DataFrame(data=np.transpose(gaussians), index=x_range, columns=names)
        self.xps.gaussians = xps_df

    def _get_binding_energies(self, specie: sym.Symbol):

        bes = []
        for orbital in self.species_manager[specie].orbitals:
            bes.append(orbital.binding_energy)
        return bes

    def _get_species_gauss(self, specie, x_range, t):

        gaussian = np.zeros(x_range.size)
        concentration = self.at(t)[specie]

        for orbital in self.species_manager[specie].orbitals:
            gaussian += concentration * orbital.splitting * \
                    stats.norm.pdf(x_range, orbital.binding_energy, _SIGMA)

        return gaussian

    # --- Plotting Methods ---------------------------------------------------

    def plot(self, exp_type='timeseries',
             species: List[sym.Symbol] = [],
             ignore: List[sym.Symbol] = [],
             t=-1, **kwargs):
        """Plot the solution as if it's a (specified) experimental observable.

        Args:
            exp_type: str, the type of experimental observable to plot as.
            species: A list of sym.Symbols. the species to plot. If None,
                will resort to what you specify with ignore.
            ignore: A list of sym.Symbols you don't want to plot. Won't matter
                if you specify species.
            t: The snapshot in time to plot (if applicable)
            **kwargs: Forwarded to DataFrame.plot/plt.plot
        """
        species = self._get_species_not_ignored(species, ignore)

        # Assume negative time means time-max
        if t < 0:
            t = self.t[-1]

        if exp_type == 'timeseries':
            self._plot_time_series(species, **kwargs)
        elif exp_type == 'spectro':
            if self.xps:
                self._plot_spectro_with_xps(species, t=t, **kwargs)
            else:
                self._plot_spectro(species, t=t, **kwargs)
        else:
            raise ValueError(f'{exp_type} not a valid plot type.')

    def _plot_time_series(self, species, **kwargs):
        """Plots the species specified as a time series."""
        self.df[species].plot(**kwargs)

    def _plot_spectro(self, species: List[sym.Symbol], t: float, **kwargs):

        binding_energies = []
        for specie in species:
            binding_energies.extend(self._get_binding_energies(specie))

        x_lower = min(binding_energies) - _PLOT_MARGIN
        x_upper = max(binding_energies) + _PLOT_MARGIN
        x_range = np.arange(x_lower, x_upper, _PLOT_RESOLUTION)

        gaussians = [self._get_species_gauss(x_range, specie, t) for specie in species]
        envelope = sum(gaussians)

        for index, gauss in enumerate(gaussians):
            name = self.species_manager[species[index]].name
            plt.fill(x_range, gauss, label=name, color=_COLORS[index])

        plt.plot(x_range, envelope, color='black', linewidth=4)

        plt.legend()
        plt.title(f'time={t}')
        plt.gca().invert_xaxis()  # Spectroscopy plots
        plt.show()

    def _plot_spectro_with_xps(self, species: List[sym.Symbol],
                               gas_interval=(), t: float = -1, **kwargs):

        # We plot on the same range that the XPS does
        x_range = self.xps.binding_energy

        gaussians = [self._get_species_gauss(x_range, specie, t) for specie in species]
        envelope = sum(gaussians)

        # Scale: max(envelope) should equal max(self.xps.intensity).
        factor = max(self.xps.intensity) / max(envelope)
        envelope *= factor
        gaussians = [gauss * factor for gauss in gaussians]

        if gas_interval:
            gas_phase = self._get_gas_phase_gauss(x_range, gas_interval)
            envelope += gas_phase

        for index, gauss in enumerate(gaussians):
            name = self.species_manager[species[index]].name
            plt.fill(x_range, gauss, label=name, color=_COLORS[index])

        if gas_interval:
            plt.fill(x_range, gas_phase, label='gas phase', color='gray')
        plt.plot(x_range, self.xps.intensity, color='green')
        plt.plot(x_range, envelope, color='black', linewidth=4)

        plt.legend()
        plt.title(f'time={t}')
        plt.gca().invert_xaxis()  # Spectroscopy plots
        plt.show()

    def _get_species_not_ignored(self, species: List[sym.Symbol] = None,
                                 ignore: List[sym.Symbol] = None):
        """A helper method to return species - ignored. If species is none,
        considers all species."""
        if not species:
            species = self.species
        if not ignore:
            ignore = []
        return [specie for specie in species if specie not in ignore]

    def _time_to_index(self, time):
        """Takes a time and returns the highest index <= that time.

        Raises:
            IndexError: If time is negative.

        Returns:
            An integer index such that self.t[index] <= time, unless no such
            index exists, in which case it returns len(self.t).
        """
        if time < 0:
            raise IndexError('time cannot be below 0.')
        return bisect.bisect_right(self.t, time) - 1

    def __repr__(self) -> str:
        shortened_sols = [sol[:10] for sol in self.states.values()]
        return 'Solution('+str(self.t[:10])+'..., '+str(shortened_sols[:10])+'...'


def simulate(rsys, time_max: float = 1, **options):
    """
    TODO
    """

    sol_t, sol_y = _solve_rsys_ode(rsys, time_max, **options)
    return Solution(sol_t, sol_y, rsys)


def _solve_rsys_ode(rsys, time_max: float = 1, **options):
    """
    TODO
    """

    ode_func = rsys.get_ode_functions()
    num_species = len(rsys._symbols)

    # schedule, a dictionary {time : [amount to add for species no. index]}
    schedule = collections.defaultdict(lambda: [0] * num_species)
    for index in range(num_species):
        for time, amount in rsys.scheduler[index].items():
            schedule[time][index] += amount

    # This is an ordered list of all the times at which we add/remove stuff.
    time_breaks = sorted(schedule.keys())

    # Do the simulation in broken pieces
    current_concs = [0] * num_species
    current_time = time_breaks.pop(0)  # Guaranteed to be 0

    while current_time < time_max:
        # Get next_time, the end of the interval we're simulating this step.
        if len(time_breaks) > 0:
            next_time = time_breaks.pop(0)
        else:
            next_time = time_max

        # Add the respective amounts for this timestep.
        for index, amount in enumerate(schedule[current_time]):
            current_concs[index] += amount

        partial_sol = integrate.solve_ivp(ode_func, (current_time, next_time), current_concs, **options)

        # Add the partial solution to the whole solution.
        if current_time == 0:
            sol_t = partial_sol.t
            sol_y = partial_sol.y
        else:
            sol_t = np.append(sol_t, partial_sol.t)
            sol_y = np.append(sol_y, partial_sol.y, axis=1)

        # Loop; set current_concs to the new ones and curren_time to the next one
        for index in range(num_species):
            current_concs[index] = sol_y[index][len(sol_y[index]) - 1]
        current_time = next_time

    # Set y for concentration equations, which the ODE solver does't calculate.
    for index, func in rsys.get_conc_functions().items():
        for tindex in range(sol_t.size):
            sol_y[index][tindex] = func(sol_t[tindex], sol_y[:, tindex])

    return sol_t, sol_y