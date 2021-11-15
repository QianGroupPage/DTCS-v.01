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

from __future__ import annotations

import bisect
from typing import Dict, List, Tuple, Union, Optional

from matplotlib import pyplot as plt
import monty.json
import numpy as np
import pandas as pd
import sympy as sym

from lblcrn.sim import bulk_crn
from lblcrn.sim import surface_crn
from lblcrn.twin import twin_abc
from lblcrn.twin import xps

from lblcrn.common import util
from lblcrn.common.colors import color_map

BulkCRNSpec: TypeAlias = 'BulkCRNSpec'


class CRNTimeSeries(twin_abc.Experiment):
    """A time series solution of an chemical reaction ODE, with utilities.

    Defines a solution to a system of differential equations, providing a
    variety of functions to manipulate and visualize the results

    Attributes:
        df: The pandas DataFrame associated with the time series.
        rsys: The reaction system which the Solution is a solution for.
        species_manager: The SpeciesManager in use.
    """
    
    def __init__(self, t: List[float], y: List[List[float]], crn):
        """Translates a t, y pair from solve_ivp into a df and initializes."""
        # TODO: Doc args
        super().__init__()

        self.df = pd.DataFrame(data=np.transpose(y),
                               index=pd.Index(t, name='time'),
                               columns=pd.Index(crn.rsys.get_symbols_ordered(),
                                                name='species'))
        self.rsys = crn.rsys
        self.species_manager = crn.species

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
        # TODO(Andrew): I don't think this style of function is transparent
        #  to the user, I'm going to force them to use xps_with().
        raise NotImplementedError()
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
        # Ignore species without binding energies
        ignore = ignore or []
        for specie in self.species:
            if specie in self.species_manager and \
                    hasattr(self.species_manager[specie], 'orbitals'):
                continue
            ignore.append(specie)

        species = twin_abc._get_species_not_ignored(species, ignore,
                                                    self.species)

        snapshot = self.at(t)
        species_concs = {}
        for specie, conc in snapshot.items():
            if specie in species:
                species_concs[specie] = conc
        if not title:
            title = f'time={snapshot.name}'
        # TODO(Andrew): Is this paradigm better than simulate_xps?
        #  I should standardize how to do it. I think that the __init__ would
        #  make more sense to have the bulk of the info and then the simulate
        #  function just forwards.
        self._xps = xps.XPSExperiment(species_manager=self.species_manager,
                                      title=title,
                                      x_range=x_range,
                                      species=species,
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
            self.df[specie].plot(ax=ax, color=color_map[specie], **kwargs)

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


class SurfaceCRNTimeSeries(CRNTimeSeries):
    """
    TODO
    It's going to have a dataframe with:
    - the time series of each run, through a multiindex
    - the map for each step in the time series (what if it's too big?)

    A property to calculate:
    - the mean result
    - the rolling mean for a window
    - the rolling mean for a run(?)
    """

    def __init__(self):
        twin_abc.Experiment.__init__(self)

    @classmethod
    def from_runs(cls, trajectories, crn):
        scts = cls()

        scts.trajs = trajectories
        scts.crn = crn
        scts.rsys = crn.rsys
        scts.species_manager = crn.species

        scts.df_full = cls.resample(scts)
        scts.df = scts.df_full['rolling_mean']
        scts._xps = None

        return scts

    def resample(self, roll_window=0.5, time_step=0.01, index=None):
        time_max = self.crn.time

        # Resample:
        sample_on = index or np.arange(0, time_max, time_step)
        time_index = set(util.flat([list(df.index) for df in self.trajs]))
        time_index.update(sample_on)

        self_df = pd.DataFrame(index=sorted(time_index))

        for index in range(len(self.trajs)):
            for col in self.trajs[index].columns:
                self_df[f'seed{index}', sym.Symbol(col)] = self.trajs[index][col]

        self_df.columns = pd.MultiIndex.from_tuples(self_df.columns)
        self_df = self_df.interpolate(method='index')
        self_df = self_df.loc[sample_on]  # TODO: Make sure to save the raw data, this drops it!

        # Calculate mean
        for species in self_df.columns.levels[1]:
            self_df['mean', species] = self_df.xs(species, axis=1, level=1).mean(axis=1)

        # Calculate rolling mean
        def rolling_mean(df, window):
            def roll(row):
                time = row.name
                return df[max(0, time - window):time].mean()
            return roll

        rolled = self_df.apply(rolling_mean(self_df['mean'], roll_window), axis=1)

        for species in rolled.columns:
            self_df['rolling_mean', species] = rolled[species]

        # Sort columns
        self_df = self_df.sort_index(axis=1)

        return self_df

    def as_dict(self):
        raise NotImplemented()

def simulate_bulk_crn(
        crn: BulkCRNSpec,
        #rsys: BulkRxnSystem = None,
        time: Optional[float] = None,
        end_when_settled: bool = False,
        #title: str = "",
        #species: List[sym.Symbol] = None,
        #ignore: List[sym.Symbol] = None,
        #x_range: Optional[np.ndarray] = None,
        #scale_factor: float = None,
        #experimental: pd.Series = None,
        #gas_interval: Tuple[float, float] = None,
        #contam_spectra: Optional[Dict[sym.Symbol, pd.Series]] = None,
        #deconv_species: Optional[List[sym.Symbol]] = None,
        #autoresample: bool = True,
        #autoscale: bool = True,
        **options
):
    """Simulate the given reaction system over time.

    Args:
        TODO: add rest
        #rsys: ReactionsSystem, the reaction system to simulate
        time: The time until which to simulate.
        species: The species to include in the XPS.
        ignore: The species to not include in the XPS.
        autoresample: Decides if the XPS resamples on edits.
        autoscale: Decides if the XPS will automatically scale the gaussians and envelope to match the experimental data.
        experimental: The experimental value of the XPS.
        gas_interval: The interval in which the peak of the gas phase is in the XPS.
        scale_factor: The scale factor by which to scale the simulated gaussians in the XPS.
        title: The name to give the XPS, used in plotting.
        options: Forwarded to scipy.integrate.solve_ivp

    Returns:
        An XPSExperiment with the simulation results as well as a CRNTimeSeries object with the time series data.
    """
    time = time or crn.time
    end_when_settled = end_when_settled or (time is None)

    sol_t, sol_y = bulk_crn.solve_rsys_ode(
        rsys=crn.rsys,
        time_max=time,
        end_when_settled=end_when_settled,
        max_step=crn.max_step,
        **options)
    cts = CRNTimeSeries(sol_t, sol_y, crn)

    return cts


def simulate_surface_crn(scrn, **kwargs):
    ens = surface_crn.scrn.scrn_simulate(scrn.rsys,
                                         time_max=scrn.time,
                                         ensemble_size=scrn.runs,

                                         video=False,
                                         spectra_in_video=False,
                                         video_path='output',
                                         trajectory_path='output',)
    results = None
    if scrn.runs == 1:
        results = [ens]
    else:
        results = ens.results

    scts = SurfaceCRNTimeSeries.from_runs([result.df_raw for result in ens.results], scrn)
    return scts
