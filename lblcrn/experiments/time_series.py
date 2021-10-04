"""Utilities for simulating reaction systems over time.

Exports:
    CRNTimeSeries: A time-series solution of a chemical reaction system.
    simulate(): Simulate a given reaction system over time.

Example:
    Once you have your reaction system rsys set up, you can proceed as
    follows:

    xps, time_series = simulate(rsys, time_max=20)
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
from typing import Dict, List, Tuple, Union

from matplotlib import pyplot as plt
import monty.json
import numpy as np
import pandas as pd
import sympy as sym

from lblcrn import bulk_crn
from lblcrn.experiments import experiment
from lblcrn.experiments import xps
from lblcrn.crn_sym.rxn_system import RxnSystem


class CRNTimeSeries(experiment.Experiment):
    """A time series solution of an chemical reaction ODE, with utilities.

    Defines a solution to a system of differential equations, providing a
    variety of functions to manipulate and visualize the results

    Attributes:
        df: The pandas DataFrame associated with the time series.
        rsys: The reaction system which the Solution is a solution for.
        species_manager: The SpeciesManager in use.
    """
    
    def __init__(self, t: List[float], y: List[List[float]], rsys):
        """Translates a t, y pair from solve_ivp into a df and initializes."""
        # TODO: Doc args
        super().__init__()

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

    def xps_with(self, t: float = -1,
                 title: str = '',
                 species: List[sym.Symbol] = None,
                 ignore: List[sym.Symbol] = None,
                 x_range: np.ndarray = None,
                 scale_factor: float = None,
                 experimental: pd.Series = None,
                 gas_interval: Tuple[float, float] = None,
                 contam_spectra: Dict[sym.Symbol, pd.Series] = None,
                 deconv_species: List[sym.Symbol] = None,
                 autoresample: bool = True,
                 autoscale: bool = True,):
        """Calculates a simulated XPS observable at time t.

        In addition to returning, saves the information into self.xps.

        Args:
            t: The time at which to take a snapshot. Defaults to the max time.
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

        Returns:
            An XPSExperiment object with the parameters you specified.
        """
        species = experiment._get_species_not_ignored(species, ignore,
                                                      self.species)
        snapshot = self.at(t)
        species_concs = {}
        for specie, conc in snapshot.items():
            if specie in species:
                species_concs[specie] = conc
        if not title:
            title = f'time={snapshot.name}'
        self._xps = xps.XPSExperiment(species_manager=self.species_manager,
                                      title=title,
                                      x_range=x_range,
                                      scale_factor=scale_factor,
                                      sim_concs=species_concs,
                                      experimental=experimental,
                                      gas_interval=gas_interval,
                                      contam_spectra=contam_spectra,
                                      deconv_species=deconv_species,
                                      autoresample=autoresample,
                                      autoscale=autoscale,)
        return self._xps

    # --- Plotting -----------------------------------------------------------

    def _plot(self, ax: plt.Axes, species: List[sym.Symbol], **kwargs):
        """Plot the reaction network time series.

        Args:
            ax: The plt.Axes on which to plot.
            species: A list of sym.Symbols, the species to plot.
            **kwargs: Forwarded.
        """

        for i, name in enumerate(species):
            if isinstance(name, str):
                species[i] = sym.Symbol(name)
            specie = species[i]
            self.df[specie].plot(ax=ax, color=self.species_manager[specie].color, **kwargs)

    # --- Utility -------------------------------------------------------------

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

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict()
        d['rsys'] = self.rsys.as_dict()
        d['t'] = self.t
        d['y'] = np.transpose(self.df.to_numpy())
        return d

    def tof(self, species: Union[str, sym.Symbol], step_difference: bool=False) -> float:
        """
        TOF series for the result

        :param species: species name in either string or sympy.Symbol
        :param step_difference: if True, define TOF as local concentration difference / local time step size;
                                otherwise, define TOF as concentration so far / time till now;
        # :param n: number used to calcualte n-th discrete difference if step_difference is True.
        :return: TOF series.
        """
        if isinstance(species, sym.Symbol):
            species = species.name
        y = [self.df[colname] for colname in self.df.columns if str(colname) == species][0]
        if len(y.index) == 0:
            print("The species has an empty time series.")
            return None
        if step_difference:
            n = 1
            tof_series = np.diff(y, n=n) / np.diff(y.index, n=n)
            time_series = y.index[:-n]
        else:
            tof_series =  y / y.index
            time_series = tof_series.index
        tof_df = pd.DataFrame.from_dict({"Time /s": time_series, "TOF /s-1": tof_series})
        tof_df.set_index("Time /s", inplace=True)
        return tof_df

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        decode = monty.json.MontyDecoder().process_decoded
        d['rsys'] = decode(d['rsys'])
        d['t'] = decode(d['t'])
        d['y'] = decode(d['y'])
        return cls(**d)
