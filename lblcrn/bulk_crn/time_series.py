"""Utilities for simulating reaction systems over time.

Exports:
    CRNTimeSeries: A time-series solution of a chemical reaction system.
    simulate(): Simulate a given reaction system over time.

Example:
    Once you have your reaction system rsys set up, you can proceed as
    follows:

    time_series = simulate_crn(rsys, time_max=20)
    time_series.df  # Shows the DataFrame of the time series.
    time_series.plot()  # Plots the network's time series.
    time_series.at(t=2)  # Shows the state near time=2

    CRNTimeSeries also has experiment-specific members. You can use them to
    show simulated observables. For example, x-ray spectroscopy:

    time_series.xps_with(ignore=[species1], t=12)
    time_series.xps.experimental = experimental_data
    time_series.xps.gas_interval = (500, 510)
    time_series.xps.plot()

    For more information, see xps.py.
"""

import bisect
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import sympy as sym

from lblcrn.bulk_crn import xps
from lblcrn.bulk_crn import common
from lblcrn.crn_sym import reaction


class CRNTimeSeries:
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
    def xps(self) -> xps.XPSExperiment:
        """Gives the xps observable you've calculating, calculating a default
        if you haven't run calc_xps yet."""
        if self._xps is None:
            self.xps_with()
        return self._xps

    def xps_with(self, t: float = -1, species: List[sym.Symbol] = None,
                 ignore: List[sym.Symbol] = None,
                 experimental: pd.Series = None,
                 gas_interval: Tuple[float, float] = None,
                 scale_factor: float = 0):
        """Calculates a simulated XPS observable at time t.

        In addition to returning, saves the information into self.xps.

        Args:
            t: The time at which to take a snapshot. Defaults to the max time.
            species: The species to include in the XPS.
            ignore: The species to not include in the XPS.
            experimental: The experimental value of the XPS.
            gas_interval: The interval in which the peak of the gas phase is
                in the XPS.
            scale_factor: The scale factor by which to scale the simulated
                gaussians in the XPS.

        Returns:
            An XPSExperiment object with the parameters you specified.
        """
        species = self._get_species_not_ignored(species, ignore)
        snapshot = self.at(t)
        species_concs = {}
        for specie, conc in snapshot.items():
            if specie in species:
                species_concs[specie] = conc
        self._xps = xps.XPSExperiment(species_concs, self.species_manager,
                                      experimental=experimental,
                                      gas_interval=gas_interval,
                                      scale_factor=scale_factor,
                                      title=f'time={snapshot.name}')
        return self._xps

    # --- Plotting -----------------------------------------------------------

    def plot(self, species: List[sym.Symbol] = [],
             ignore: List[sym.Symbol] = [], **kwargs):
        """Plot the reaction network time series.

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

    def _repr_html_(self) -> Optional[str]:
        """For iPython, mostly just gives the DataFrame."""
        # TODO: do this without accessing private member?\
        df_html = self.df.head()._repr_html_()

        html = f"""<pre>{self.__class__.__name__} with head:</pre>
        {df_html}"""

        return html


def simulate_crn(rsys: reaction.RxnSystem, time_max: float = 1,
                 **options) -> CRNTimeSeries:
    """Simulate the given reaction system over time.

    Args:
        rsys: ReactionsSystem, the reaction system to simulate
        time_max: The time until which to simulate.
        **options: Forwarded to scipy.integrate.solve_ivp

    Returns:
        A CRNTimeSeries object with the concentrations over time.
    """

    sol_t, sol_y = common.solve_rsys_ode(rsys, time_max, **options)
    return CRNTimeSeries(sol_t, sol_y, rsys)


