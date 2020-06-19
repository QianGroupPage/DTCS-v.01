"""Utilities for simulating reaction systems over time.

Exports:
    Solution: A time-series solution of a chemical reaction system.
    simulate(): Simulate a given reaction system over time.

Example:
    Once you have your reaction system rsys set up, you can proceed as
    follows:

    solution = simulate(rsys, time_max=20)
    solution.df  # Shows the DataFrame of the solution.
    solution.plot()  # Plots the solution as a time series.
    solution.at(t=2)  # Shows the state near time=2

    Solution also has experiment-specific members. You can use them to show
    simulated observables. For example, x-ray spectroscopy:

    solution.calc_xps(ignore=[species1], t=12)
    solution.xps.experimental = experimental_data
    solution.xps.gas_interval = (500, 510)
    solution.xps.plot()

    For more information, see xps.py.
"""

import bisect
import collections
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy import integrate
import sympy as sym

from lblcrn.bulk_crn import xps


class Solution:
    """A time series solution of an chemical reaction ODE, with utilities.

    Defines a solution to a system of differential equations, providing a
    variety of functions to manipulate and visualize the results

    Attributes:
        df: The pandas DataFrame associated with the time series.
        rsys: The reaction system which the Solution is a solution for.
        species_manager: The SpeciesManager in use.
    """
    
    def __init__(self, t: List[float], y: List[List[float]], rsys):

        self.df = pd.DataFrame(data=np.transpose(y),
                               index=pd.Index(t, name='time'),
                               columns=pd.Index(rsys.get_symbols(),
                                                name='species'))
        self.rsys = rsys
        self.species_manager = rsys.species_manager

        self._xps = None

    # --- Accessors ----------------------------------------------------------

    def at(self, t: float):
        """Gives the state of the system at time <= t, or at end if t < 0."""
        if t < 0:
            t = self.t[-1]
        return self.df.iloc[self._time_to_index(t)]

    @property
    def species(self) -> List[sym.Symbol]:
        """Give the species' symbols in this solution."""
        return self.rsys.get_symbols()

    @property
    def t(self) -> np.ndarray:
        """Give the time axis of the timeseries."""
        return np.asarray(self.df.index)

    # --- Experiment Simulation -----------------------------------------------

    @property
    def xps(self) -> xps.XPSObservable:
        """Gives the xps observable you've calculating, calculating a default
        if you haven't run calc_xps yet."""
        if self._xps is None:
            self.calc_xps()
        return self._xps

    def calc_xps(self, t: float = -1, species: List[sym.Symbol] = None,
                 ignore: List[sym.Symbol] = None,
                 experimental: pd.Series = None,
                 gas_interval: Tuple[float, float] = None,
                 scale_factor: float = 0,
                 ):
        """Calculates a simulated XPS observable at time t.

        Doesn't return; instead, saves the information into self.xps.

        Args:
            t: The time at which to take a snapshot. Defaults to the max time.
            species: The species to include in the XPS.
            ignore: The species to not include in the XPS.
            experimental: The experimental value of the XPS.
            gas_interval: The interval in which the peak of the gas phase is
                in the XPS.
            scale_factor: The scale factor by which to scale the simulated
                gaussians in the XPS.
        """
        species = self._get_species_not_ignored(species, ignore)
        snapshot = self.at(t)
        species_concs = {}
        for specie, conc in snapshot.items():
            if specie in species:
                species_concs[specie] = conc
        self._xps = xps.XPSObservable(species_concs, self.species_manager,
                                      experimental=experimental,
                                      gas_interval=gas_interval,
                                      scale_factor=scale_factor,
                                      title=f'time={snapshot.name}')

    # --- Plotting -----------------------------------------------------------

    def plot(self, species: List[sym.Symbol] = [],
             ignore: List[sym.Symbol] = [], **kwargs):
        """Plot the solution as a time-series.

        Args:
            species: A list of sym.Symbols. the species to plot. If None,
                will resort to what you specify with ignore.
            ignore: A list of sym.Symbols you don't want to plot. Won't matter
                if you specify species.
            **kwargs: Forwarded to DataFrame.plot/plt.plot
        """
        species = self._get_species_not_ignored(species, ignore)
        self.df[species].plot(**kwargs)

    # --- Utility -------------------------------------------------------------

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

    def _repr_html(self):
        """For iPython, mostly just gives the DataFrame."""
        # TODO: do this without accessing private member?
        # TODO: Also it doesn't appera to get called.
        return self.df._repr_html_()


def simulate(rsys, time_max: float = 1, **options):
    """Simulate the given reaction system over time.

    Args:
        rsys: ReactionsSystem, the reaction system to simulate
        time_max: The time until which to simulate.
        **options: Forwarded to scipy.integrate.solve_ivp

    Returns:
        A Solution object describing the solution.
    """

    sol_t, sol_y = _solve_rsys_ode(rsys, time_max, **options)
    return Solution(sol_t, sol_y, rsys)


def _solve_rsys_ode(rsys, time_max: float = 1, **options):
    """Simulate the given reaction system over time.

    Private in case we want to add multiple possible solving methods.

    Args:
        rsys: ReactionsSystem, the reaction system to simulate
        time_max: The time until which to simulate.
        **options: Forwarded to scipy.integrate.solve_ivp

    Returns:
        A Solution object describing the solution.
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

        partial_sol = integrate.solve_ivp(ode_func, (current_time, next_time),
                                          current_concs, **options)

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