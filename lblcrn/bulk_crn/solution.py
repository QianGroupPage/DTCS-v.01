import bisect
import collections
from typing import List, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import integrate
from scipy import stats
from sklearn import metrics
import sympy as sym

class Solution:
    """A time series solution of an chemical reaction ODE, with utilities.

    Defines a solution to a system of differential equations, providing a
    variety of functions to manipulate and visualize the results

    Attributes:
        df: The pandas DataFrame associated with the time series.
    """
    
    def __init__(self, t: List[float], y: List[List[float]], rsys):

        self.df = pd.DataFrame(data=np.transpose(y),
                               index=pd.Index(t, name='time'),
                               columns=pd.Index(rsys.get_symbols(),
                                                name='species'))
        self.rsys = rsys
        self.species_manager = rsys.species_manager
        self._default_ignore = []

        # TODO ---------------------------------------------------------------

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
        self.xps = None

    # --- Accessor Methods ---------------------------------------------------

    def at(self, t: float):
        """Gives the state of the system at time <= t, or at end if t < 0."""
        if t < 0:
            t = self.t[-1]
        return self.df.iloc[self._time_to_index(t)]
    
    def at_ini(self):
        """Gives the initial state."""
        return self.at(0)
    
    def at_end(self):
        """Gives the end state."""
        return self.at(-1)

    @property
    def species(self) -> List[sym.Symbol]:
        """Give the species' symbols in this solution."""
        return self.rsys.get_symbols()

    # --- Plotting Methods ---------------------------------------------------

    def ignore(self, species: List[sym.Symbol]):
        """Mark species to be ignored by the plotter by default."""
        self._default_ignore = species

    def _not_ignored(self, ignored: List[sym.Symbol]) -> List[sym.Symbol]:
        """Get all the species not in ignored"""
        return [symbol for symbol in self.species if symbol not in ignored]

    def plot(self, exp_type='timeseries',
             species: List[sym.Symbol] = None,
             ignore: List[sym.Symbol] = None,
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
        # Decide which species to plot
        if species:
            pass
        elif ignore:
            species = self._not_ignored(ignore)
        else:
            species = self._not_ignored(self._default_ignore)

        # Assume negative time means time-max
        if t < 0:
            t = self.t[-1]

        if exp_type == 'timeseries':
            self._plot_time_series(species, **kwargs)
        elif exp_type == 'spectro':
            if self.xps:
                self._plot_spectro_with_xps(species, t, **kwargs)
            else:
                self._plot_spectro(species, t, **kwargs)
        else:
            raise ValueError(f'{exp_type} not a valid plot type.')

    def _plot_time_series(self, species, **kwargs):
        """Plots the species specified as a time series."""
        self.df[species].plot(**kwargs)

    # --- Spectro Plotting Methods --------------------------------------------

    def _plot_spectro(self, species: List[sym.Symbol], t: float, **kwargs):
        """Plots spectrum without XPS

        Args:
            species:
            t:
            **kwargs:
        """
        MARGIN = 5
        RESOLUTION = 0.001
        COLORS = ['red', 'green', 'orange', 'blue', 'purple', 'pink', 'yellow', 'gray', 'cyan']

        # Gets binding energy data, format {species: [(energy, splitting), ...]}
        binding_energies = self._get_binding_energies(species)

        all_energies = []
        for orbitals in binding_energies.values():
            for orbital in orbitals:  # These are not orbital objects, just pairs.
                all_energies.append(orbital[0])

        x_lower = min(all_energies) - MARGIN
        x_upper = max(all_energies) + MARGIN
        x_range = np.arange(x_lower, x_upper, RESOLUTION)

        # Get that what I am going to plot.
        gaussians, envelope = self._calculate_gaussians(species, binding_energies, x_range, t)

        # Plot all the gaussians
        for index, gauss in enumerate(gaussians):
            name = self.species_manager[species[index]].name
            plt.fill(x_range, gauss, label=name, color=COLORS[index])

        # Plot the envelope
        plt.plot(self.resampled_binding_energies, self.envelope, linewidth=4, color='black')

        plt.legend()
        plt.gca().invert_xaxis()
        plt.show()

    def _plot_spectro_with_xps(self, species: List[sym.Symbol], t: float, **kwargs):
        pass

    def _get_binding_energies(self, species: List[sym.Symbol], bounds=()):
        """Get the binding energies for the given species.

        Args:
            species: A list of sym.Symbols, the species to get the binding
                energies from.
            bounds: Optionally restrict energies to interval (lower, upper).

        Returns:
            A dictionary {symbol: binding energy info} for each specie given.
            The binding energy info is a list of pairs (energy, splitting).
        """
        bes = collections.defaultdict(list)
        for symbol in species:
            for orbital in self.species_manager[symbol].orbitals:
                if not bounds or bounds[0] < orbital.binding_energy < bounds[1]:
                    bes[symbol].append((orbital.binding_energy, orbital.splitting))
        return bes

    def _calculate_gaussians(self, species, binding_energies, x_range, t):

        SIGMA = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)

        # The envelope is the sum of all gaussians
        envelope = np.zeros(x_range.size)
        # A list of gaussians
        gaussians = []
        concentrations = self.at(t)

        print(binding_energies)

        for specie in species:
            for energy, splitting in binding_energies[specie]:
                gauss = concentrations[specie] * splitting * stats.norm.pdf(x_range, energy, SIGMA)
                envelope += gauss
                gaussians.append(gauss)

        return gaussians, envelope

# TODO ------------------------------------------------------------------- (TEMP) --------------------------------

    def plot_gaussian(self, envelope: bool = False, overlay: bool = False, resample_envelope: bool =
    False, ax=None, title=''):
        """
        Plots a gaussian distribution of the final species concentrations. FWHM is set at 0.75
        If specified, an envelope curve is also plotted
        """
        colors = ['red', 'green', 'orange', 'blue', 'purple', 'pink', 'yellow', 'gray', 'cyan']
        if not ax:
            for i, dist in sorted(enumerate(self.distributions), key=lambda x: max(x[1]), reverse=True):
                plt.fill(self.binding_energies, dist, label=self.names[i], color=colors[i])
            plt.legend()

            if overlay:
                if resample_envelope:
                    plt.plot(self.resampled_binding_energies, self.resampled_intensity, color='green')
                else:
                    plt.plot(self.xps.binding_energy, self.xps.intensity, color='green')

            if envelope:
                plt.plot(self.resampled_binding_energies, self.envelope, linewidth=4, color='black')

            plt.gca().invert_xaxis()
            plt.show()
        else:
            for i, dist in sorted(enumerate(self.distributions), key=lambda x: max(x[1]), reverse=True):
                ax.fill(self.binding_energies, dist, label=self.names[i], color=colors[i])
            ax.legend()

            if overlay:
                if resample_envelope:
                    ax.plot(self.resampled_binding_energies, self.resampled_intensity, color='green')
                else:
                    ax.plot(self.xps.binding_energy, self.xps.intensity, color='green')

            if envelope:
                ax.plot(self.resampled_binding_energies, self.envelope, linewidth=4, color='black')

            ax.set_xlim(max(self.resampled_binding_energies), min(self.resampled_binding_energies))
            ax.title.set_text(title)

    def scale(self, to_scale, exp):
        """Scale experimental data intensity to match that of the experimental data.

        The largest value of the simulated data is scaled to match the largest value of the
        experimental data, and this scaling factor is then applied across all simulated data.
        """
        max_to_scale = max(to_scale[1])
        max_exp = max(exp[1])
        new_envelope = []
        new_dists = []
        scaling = max_exp / max_to_scale

        for v in to_scale[1]:
            new_envelope.append(v * scaling)

        for d in to_scale[2]:
            new_dists.append(d * scaling)

        return new_envelope, new_dists

    def resample(self):
        """Resample the simulated data, reducing its size to match that of the experimental data.
        """
        # rei = []
        # i = 0
        # bes = self.binding_energies
        # for intensity, be in zip(list(reversed(self.xps.intensity)), list(reversed(self.xps.binding_energy))):
        # while i < len(bes) and be > bes[i]:
        # rei.append(intensity)
        # i += 1
        # if i >= len(bes):
        # break
        # self.resampled_intensity = rei
        e = self.envelope
        bes = self.binding_energies
        xbes = self.xps.binding_energy
        xi = self.xps.intensity

        r_e, r_bes = [], []
        i = 0
        for b in list(reversed(xbes)): # TODO: reverse it automatically
            while i < len(e) - 1 and bes[i] < b:
                i += 1
            r_e.append(e[i])
            r_bes.append(b)
        self.envelope = r_e
        self.resampled_binding_energies = np.array(r_bes)
        self.resampled_intensity = list(reversed(xi))

    def process(self, gas_range=()):
        """Resample, scale, and calculate envelopes and other characteristics of the data.

        First the binding energy bounds are found, the envelope curve and individual species
        gaussians are then computed. Finally, if an experimental data object has been set, the
        simulated data is resampled and sacled.
        """
        SIGMA = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)

        min_be = float('inf')
        max_be = float('-inf')

        # determine x axis bounds
        for name, substance in self.substances.items():
            if name not in self._default_ignore:
                min_be = min(substance.orbitals[0].binding_energy, min_be)
                max_be = max(substance.orbitals[0].binding_energy, max_be)

        self.binding_energies = np.arange(min_be - 5, max_be + 5, .001)
        self.envelope = np.zeros(self.binding_energies.size)
        self.distributions = []
        self.names = []

        for name, sol in self.final_state().items():
            if name not in self._default_ignore:
                for o in self.substances[name].orbitals:
                    be = o.binding_energy
                    dist = sol * stats.norm.pdf(self.binding_energies, be, SIGMA)
                    self.envelope += dist
                    self.distributions.append(dist)
                    self.names.append(name)

        self.resampled_binding_energies = self.binding_energies

        if self.xps:
            self.resample()
            self.envelope, self.distributions = self.scale((self.resampled_binding_energies, self.envelope,
                                                            self.distributions),
                                                           (self.xps.binding_energy, self.xps.intensity))

            # If a valid gas range is specified, create a fake peak
            if len(gas_range) == 2:
                start = gas_range[0]
                end = gas_range[1]
                intensities = list(reversed(self.xps.intensity))
                bes = list(reversed(self.xps.binding_energy))
                print(bes)

                i = 0
                while i < len(bes) and bes[i] < start:
                    i += 1

                if i < len(bes):
                    max_intensity = intensities[i]
                    be = i
                    while i < len(bes) and i < bes[i] < end:
                        if intensities[i] > max_intensity:
                            max_intensity = intensities[i]
                            be = bes[i]
                        i += 1
                    if i < len(bes):
                        dist = max_intensity * stats.norm.pdf(self.binding_energies, be, SIGMA)
                        resampled_dist = max_intensity * stats.norm.pdf(self.resampled_binding_energies, be, SIGMA)
                        self.envelope += resampled_dist
                        self.distributions.append(dist)
                        self.names.append('gas phase')

# TODO ------------------------------------------------------------------- (TEMP) --------------------------------

    def set_experimental(self, xps):
        """Takes an xps object contain experimental data, and scales the simulated solution.
        """
        self.xps = xps

    def sols(self) -> List[List[float]]:  # TODO
        return list(self.states.values())

    def time_steps(self) -> List[float]:  # TODO
        return self.t
    
    def basic_plot(self):  # TODO
        """Draws a basic plot over the time domain
        """
        for name, sol in self.states.items():
            plt.plot(self.t, sol, label=name)
        plt.legend()

    def final_state(self) -> Dict[str, float]: # TODO
        """Returns the final state of the system of equations
        """
        final_state: Dict[str, float] = {}
        for name, sol in self.states.items():
            final_state[name] = sol[len(sol) - 1]
        return final_state

    def var_sols(self, *vars) -> Dict[str, List[float]]:  # TODO
        """Returns a list of solutions for the specified variables passed in as arguments
        """
        return {name: sol for name, sol in self.states.items() if name in vars}

    def rmse(self):
        return np.sqrt(metrics.mean_squared_error(self.resampled_intensity, self.envelope))

    def mae(self):
        return metrics.mean_absolute_error(self.resampled_intensity, self.envelope)
    
    def integral_diff(self):
        return abs(integrate.trapz(self.resampled_intensity, self.resampled_binding_energies) -
                integrate.trapz(self.envelope, self.resampled_binding_energies))

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

    return sol_t, sol_y