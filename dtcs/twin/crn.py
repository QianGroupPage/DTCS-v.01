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

from typing import Dict, List, Tuple, Union, Optional, Collection

import bisect
import os

from matplotlib import pyplot as plt
from matplotlib import gridspec
import monty.json
import numpy as np
import pandas as pd
import sympy as sym

from dtcs.io import scrn_video
from dtcs.twin import twin_abc
from dtcs.twin import xps


from dtcs.common import util

from dtcs.common.display import color_map

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
        self.species_manager = crn.sm

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

    @property
    def time_max(self) -> float:
        """Give the duration of the simulation"""
        return self.crn.time

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

    def _plot(self,
              ax: plt.Axes,
              species: List[sym.Symbol],
              legend: bool = True,
              t_lines: Optional[Union[float, Tuple]] = None,
              **kwargs):
        """Plot the reaction network time series.

        Args:
            ax: The plt.Axes on which to plot.
            species: A list of sym.Symbols, the species to plot.
            **kwargs: Forwarded.
        """

        # Plot each species
        for i, name in enumerate(species):
            if isinstance(name, str):
                species[i] = sym.Symbol(name)
            specie = species[i]
            self.df[specie].plot(ax=ax, color=color_map[specie], **kwargs)

        if legend:
            ax.legend()

        if t_lines:
            if not isinstance(t_lines, Collection):
                t_lines = (t_lines, )
            for time in t_lines:
                ax.axvline(x=time, color='black')

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
        scts.species_manager = crn.sm

        scts.df_full = cls.resample(scts)
        scts.df = scts.df_full['rolling_mean']
        scts._xps = None

        return scts

    @classmethod
    def from_times_and_counts(cls, times, species_counts, crn):
        dfs = [pd.DataFrame(species_counts[index], times[index])
               for index in range(len(times))]

        return cls.from_runs(dfs, crn)

    def at(self, t: float, run: Optional[int] = None):
        """Gives the state of the system at time <= t, or at end if t < 0.

        If you supply a run number, it will give you information from just
        that run. Otherwise, it will give you the averaged data.
        """
        if run is not None:
            df = self.df_full[f'seed{run}']
        else:
            df = self.df

        if t < 0:
            t = self.t[-1]
        return df.iloc[self._time_to_index(t)]

    def xps_with(self, t: float = -1, run: Optional[int] = None,
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
            run: Which run you want to get an XPS of, or None, for averaged
                data.
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

        snapshot = self.at(t, run=run)
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

    @staticmethod
    def _plot_video(
            surf_img, scts, run, time,
            figsize=(12, 5),
            legend=True,
            title=None,
            spectrum_kwargs=None,
            timeseries_kwargs=None,
    ):
        # --- Default Arguments ---
        spectrum_kwargs = dict(
            legend=False,
            title=False,
            peak_lines=False,
        ) | (spectrum_kwargs or {})
        timeseries_kwargs = dict(
            legend=False
        ) | (timeseries_kwargs or {})

        # --- Create a figure and a complicated gridspec for subplots ---
        # Calculate aspect ratio to not distort image
        surf_height, surf_width, _ = surf_img.shape
        fig_width, fig_height = figsize
        width_ratio = max(1 / (surf_height / surf_width * fig_width / fig_height - 1), 1)

        # Make the figure and the grid specification
        fig = plt.figure(figsize=figsize, tight_layout=True)
        gs = gridspec.GridSpec(2, 2, width_ratios=[width_ratio, 1])

        # Add the subplots
        ax_surf = fig.add_subplot(gs[:, 0])
        ax_spectrum = fig.add_subplot(gs[0, 1])
        ax_timeseries = fig.add_subplot(gs[1, 1])

        # --- Display the image of the surface ---
        ax_surf.imshow(surf_img)
        ax_surf.axis('off')

        # --- Display the XPS Spectrum ---
        scts.xps_with(
            run=run,
            t=time,
        ).plot(ax=ax_spectrum, **spectrum_kwargs)
        ax_spectrum.spines['top'].set_visible(False)
        ax_spectrum.spines['right'].set_visible(False)

        # --- Display the time series ---
        scts.plot(ax=ax_timeseries, t_lines=(time,), **timeseries_kwargs)
        ax_timeseries.spines['top'].set_visible(False)
        ax_timeseries.spines['right'].set_visible(False)

        # --- Title and legend and so on ---
        if title is None: title = 't={0:.2f}'
        if title:
            fig.suptitle(title.format(time), fontsize=16)

        if legend:
            legend_patches = util.get_legend_patches(scts.species)
            fig.legend(handles=legend_patches, loc='upper left')

        # --- Return ---
        return fig, [ax_surf, ax_spectrum, ax_timeseries]

    @util.feature('scrn-image')
    def plot_image(
            self,
            run: int = 0,
            time: int = -1,
            plot: Optional[callable] = None,
            surface_img_dpi: float = 200,
            **plot_kwargs,
    ):
        # scrn_video.make_scrn_images doesn't handle t=-1
        time = time if time > 0 else self.time_max

        # Default plotting function
        if not plot: plot = self._plot_video

        images, width, height = scrn_video.make_scrn_images(
            scts=self,
            run=run,
            frame_times=(time,),
            dpi=surface_img_dpi,
        )

        return plot(
            surf_img=images[time],
            scts=self,
            run=run,
            time=time,
            **plot_kwargs,
        )

    @util.feature('scrn-video')
    def make_video(
            self,
            run: int = 0,

            frames_dir: str = 'frames',  # TODO(Andrew) option to delete when done?
            output_fname: str = 'output',
            frames_per_timestep=10,
            frames_per_second=None,
            surface_img_dpi=200,
            fourcc='vp09',

            plot_func: Optional[callable] = None,
            **plot_kwargs
    ):
        # Default plotting function
        if not plot_func: plot_func = self._plot_video
        # Default to realtime
        if frames_per_second is None: frames_per_second = frames_per_timestep

        output_path = scrn_video.make_scrn_video(
            scts=self,
            plot_func=plot_func,
            run=run,
            frames_dir=frames_dir,
            output_fname=output_fname,
            frames_per_timestep=frames_per_timestep,
            frames_per_second=frames_per_second,
            surface_img_dpi=surface_img_dpi,
            fourcc=fourcc,
            **plot_kwargs,
        )

        print(f'Wrote to {os.path.relpath(output_path)}')
        try:
            from IPython import display

            return display.Video(os.path.relpath(output_path))
        except ModuleNotFoundError:
            return os.path.relpath(output_path)

    def as_dict(self):
        raise NotImplementedError()
