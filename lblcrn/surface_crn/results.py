import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.stats import norm
import lblcrn.surface_crn.xps as xps
import os
from IPython.display import HTML
import ffmpeg
from lblcrn.common import color_to_HEX
import matplotlib.backends.backend_agg as agg
# Set default font
matplotlib.rcParams["font.family"] = "arial"


class Results:
    """
    Defines a solution to a system of differential equations, providing a variety of
    functions to manipulate and visualize the results
    """
    TIME_COL_NAME = "Time (s) "  # column name used for the time entry

    def __init__(self, manifest_file, rxns, df=None):
        """
        :param manifest_file:
        :param rxns:
        :param df: a df whose index may be unnamed but will be reset into TIME_COL_NAME
        """

        self.manifest_file = manifest_file
        self.rxns = rxns

        self.df = df.copy()
        self.df.rename_axis(Results.TIME_COL_NAME, inplace=True)

        # Correct the count values for larger species
        for s in rxns.species_manager.large_species:
            if s.name in self.df.columns:
                # print(f"dividing species {s.name} by {s.size}")
                self.df[s.name] = self.df[s.name] / s.size

        self.resample_evolution()

        self.xps_df = None  # The xps dataframe we have for reference

        self.species_ordering = None  # The ordering of species as they appear in our figures.
        self.species_colors = {}  # A dictionary from species names to their colors

        color_index = rxns.get_colors()
        sub_s_list = []
        [sub_s_list.extend(l) for l in rxns.species_manager.to_sum_dict.values()]
        primary_s = rxns.species_manager.to_sum_dict.keys()
        species_tracked = sorted(list(set(list(rxns.get_symbols()) + list(primary_s))), key=lambda s: str(s))
        self.species_ordering = [s for s in species_tracked if s in primary_s or s not in sub_s_list]
        self.species_colors = {s: color_index[s] for s in self.species_ordering}
        self.species_ordering = [str(s) for s in self.species_ordering]
        self.species_colors = {str(s): ''.join(color_to_HEX(c)) for s, c in self.species_colors.items()}
        import sympy as sym
        self.substances = {s: rxns.species_manager.species_from_symbol(sym.Symbol(s)) for s in self.species_ordering}\
            if rxns else {}
        # species_tracked = list(rxns.get_symbols() + list(primary_s))
        # self.species_ordering = [rxns.species_manager.species_from_symbol(s) for s in species_tracked
        #                          if rxns.species_manager.symbol_in_sm(s) and (s in primary_s or s not in sub_s_list)]
        # self.species_colors = {s: color_index[s] for s in self.species_ordering}
        # self.substances = {s.name: s for s in self.species_ordering}
        # self.species_ordering = [s.name for s in self.species_ordering]
        # self.species_colors = {s.name: color_to_HEX(c) for s, c in self.species_colors.items()}

        sub_s = rxns.species_manager.to_sum_dict
        self.sum_sub_species(sub_s)
        self.sum_same_be()


        # Play the videos
        self.video = None
        self.video_trajectory = None


    # def extend(self, df):
    #     """
    #     A
    #     :param df:
    #     :return:
    #     """


    @staticmethod
    def side_by_side_axes(num_axes=2, output_fig=False):
        """
        A function designed to wrap around matplotlib to encourage
        people to do side by side comparisons more often.
        :param num_axes: number of axes to generate
        :param output_fig: output fig, axes if set to True
        :return: a list of axes
        """
        num_lines = num_axes // 2 + num_axes % 2
        fig, axes = plt.subplots(num_lines, 2, figsize=(16, 6 * num_lines))

        def trim_axs(axs, N):
            """little helper to massage the axs list to have correct length..."""
            """https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/markevery_demo.html#sphx-glr-gallery-lines-bars-and-markers-markevery-demo-py"""
            axs = axs.flat
            for ax in axs[N:]:
                ax.remove()
            return axs[:N]

        if output_fig:
            return fig, trim_axs(axes, num_axes)
        else:
            return trim_axs(axes, num_axes)

    def play_video(self, slowdown_factor=1):
        """
        Play the simulation video if a video is produced.


        :slowdown_factor: 1 by default, 30 is suggested for reasonable viewing by humans.
        :return: HTML object for playing a video in the IPython Notebook
        """
        if self.video is None:
            # TODO: does absolute path work?
            # self.video = "/Users/ye/Desktop/lbl-crn/Surface CRN Videos/scrn simulation.mp4"
            # self.video = "Surface CRN Videos/scrn simulation.mp4"
            raise Exception("There is no video generated for the reaction system.")

        if slowdown_factor == 1:
            video = self.video
        else:
            # Use a factor for slowing down the video
            head, name = os.path.split(self.video)
            name, ext = os.path.splitext(name)
            slowmo_name = f'{head}/{name}_{slowdown_factor}x_slower{ext}'
            # Overwrite existing files.
            if os.path.isfile(slowmo_name):
                os.remove(slowmo_name)
            ffmpeg.input(self.video).setpts(f"{slowdown_factor}*PTS").output(
                f'{head}/{name}_{slowdown_factor}x_slower{ext}').run()
            video = slowmo_name
        return HTML(f"""
        <video width="960" height="720" controls>
          <source src="{video}" type="video/mp4">
        </video>
        """)

    def resample_evolution(self, round=1):
        df = self.df.copy()
        df["Time (s)"] = df.index
        df["Time (s)"] = df["Time (s)"].round(round)
        df = df.groupby("Time (s)").mean()
        self.df_raw = self.df
        self.df = df.reset_index().set_index("Time (s)")

    @staticmethod
    def from_directory(path, rxns):
        """
        Construct a results object from based on a directory.
        :param path:
        :param rxns:
        :return:
        """
        if os.path.isfile(f"{path}/Data.gz"):
            return Results(f"{path}/reaction_rules.txt", rxns, pd.read_csv(f"{path}/Data.gz"))
        else:
            return Results(f"{path}/reaction_rules.txt", rxns, pd.read_csv(f"{path}/Data.csv"))

    @staticmethod
    def from_concs_times(manifest_file, rxns, concs, times):
        """
        :param manifest_file: The manifest file corresponding to the concs dictionary;
        :param concs: a dictionary from species name (string) to a list of concentration values;
        :param times: a list of time steps corresponding to the time steps for each list in concs
        :return: a new Solution object
        """
        r = Results(manifest_file, rxns, Results.concs_times_df(concs, times))
        return r

    @staticmethod
    def from_counts(rxns, counter):
        """
        :param rxns: a rxn_system object
        :param counter: a dictionary from species name to counts
        :return: a results object representing only 1 timestep.
        """
        concs = {s: [c] for s, c in counter.items()}
        for s in rxns.get_symbols():
            if str(s) not in concs:
                concs[str(s)] = [0]
        return Results.from_concs_times(None, rxns, concs, [0])

    @staticmethod
    def concs_times_df(concentrations, times):
        """
        :param concentrations: a dictionary from species name (string) to a list of concentration values;
        :param times: a list of time steps corresponding to the time steps for each list in concs
        :return: a pandas DataFrame where each
        """
        if not concentrations:
            raise Exception("The concentrations dictionary cannot be empty.")
        else:
            for key, value in concentrations.items():
                if len(value) != len(times):
                    raise Exception(f"{key} entry in the concentrations dictionary" +
                                    f" has {len(value)} number of entries whereas times list has " +
                                    f"{len(times)} number of entries.")
        df = pd.DataFrame.from_dict(concentrations)
        # TODO: get rid of "(s) ", a typo that persists from old simulations.
        df.insert(0, Results.TIME_COL_NAME, times)
        df = df.set_index(Results.TIME_COL_NAME)
        return df

    def compute_converged_values(self, convergence_time):
        """
        Average all the values after time step convergence_time.
        """
        return self.average_values(convergence_time, self.df.index[-1])

    def average_values(self, start, end=-1):
        """
        Average the result df from time start to time end.
        :return: a df with time column dropped
        """
        if end == -1:
            end = self.df_raw.index[-1]
        df = self.df_raw

        # print("Entire df", df)
        # # TODO: use the raw df
        # print("Begins the df chosen for averaging")
        #
        # print(f"\n start {start}, end {end}")
        #
        # print(df.loc[(start <= df.index) & (df.index <= end)])
        # print("Ends the df chosen for averaging")
        return df.loc[(start <= df.index) & (df.index <= end)].mean()

    def concentration_bar_plot(self, t=-1, scale=None):
        """
        Plot a concentration at time t as a bar plot. By default, plot the last time frame
        :return: None
        """
        if t == -1:
            t = self.df[self.TIME_COL_NAME].iloc[-1]
        avg_series = self.average_values(t-0.5, t+0.5)
        ax = sns.barplot(x=avg_series.index, y=avg_series)

        if scale == "Entire DF":
            ax.set_ylim(0, self.df.to_numpy().max())

        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

    def rank_by_final_value(self):
        """
        Return a list of species ranked by their final concentration value.

        Final concentration, for the purposes here, is defined as the average of the last 100 values in a time series.
        """
        species = list(self.df.columns)
        # TODO: automatically find out the convergence time.
        species.sort(reverse=True, key=lambda k: sum(self.df[k][-100:]) / 100)
        return species

    # TODO: test this.
    def sum_sub_species(self, sub_species_dict):
        """
        Produce the total concentration for each species in sub_species_dict.

        Append the resulting series of total concentration to concs.
        """
        self.df = self.df_raw
        for k, v in sub_species_dict.items():
            # Transform from symbol to strs
            k = str(k)

            v = [str(s) for s in v]

            summed_series = pd.Series(0, self.df.index)
            for s in v:
                # If parent is in the list of children, but not in the data column,
                # assume that parent is not a valid species in the underlying system.
                # TODO: improve the logic here
                if s not in self.df.columns:
                    continue
                if s not in self.df.columns and s == k:
                    continue
                if s != k and s not in self.df.columns:
                    raise Exception(f"Subspecies {s} has not been recorded in the results data frame.")
                if s in self.df.columns:
                    summed_series += self.df[s]
            self.df[k] = summed_series
        self.resample_evolution()

    def sum_same_be(self):
        """
        Sum the number of species with the same binding energy.
        :return:
        """
        return
        bes_to_names = {}
        to_sum = {}
        for name in self.species_ordering:
            for o in self.substances[name].orbitals:
                be = o.binding_energy
                if be in to_sum:
                    to_sum[be].append(name)
                elif be in bes_to_names:
                    to_sum[be] = bes_to_names[name] + [name]
                else:
                    bes_to_names[be] = [name]
        # for v in to_sum.values:
        # TODO: how shall we name the summed species?

    # TODO: decrease the figure size in case zoom = False
    def plot_evolution(self, species_in_figure=None, start_time=0, end_time=-1, title="", ax=None, save=False,
                       return_fig=False, path="", use_raw_data=False, zoom=False, show_fig=True, x_axis_xlim=0):
        """
        Plot the concentrations from start_time until time step end_time. -1 means till the end.

        Name is the specific name of this plot.
        """
        plt.style.use('seaborn-white')
        if end_time == -1:
            end_time = self.df_raw.index.max()

        # print("end time", '{:.15f}'.format(end_time))

        if use_raw_data:
            df = self.df_raw
        else:
            df = self.df
        df = df[(df.index <= end_time) & (df.index >= start_time)]

        if species_in_figure is None:
            species_in_figure = self.species_ordering

        if zoom:
            irange = range(len(species_in_figure) - 1)
        else:
            irange = [0]

        if not ax:
            ax_given = False
            fig, axes = plt.subplots(len(irange), 1)
            # TODO: fix this for multiple axes
            fig.set_figheight(6)
            fig.set_figwidth(8)
            # fig.subplots_adjust(top=0.95)
            # fig.suptitle(f"{title}", fontsize=48)
        else:
            ax_given = True
            axes = [ax]

        if not isinstance(axes, list):
            axes = [axes]

        for i in irange:
            ax = axes[i]
            ax.set_xlim(left=x_axis_xlim, right=end_time)
            # print(f"left {0}, right {end_time*1.1}")
            for j in range(i, len(species_in_figure)):
                species = species_in_figure[j]
                ax.tick_params(axis='both', which='both', labelsize=12)
                ax.plot(df[species], color=self.species_colors[species], label=species, linewidth=2)
                ax.legend(fontsize=12, numpoints=30)

            ax.set_title(title)
            ax.set_xlabel("Time (s)", fontsize=12)
            ax.set_ylabel("Molecule Count (#)", fontsize=12)

        if not show_fig:
            if ax_given:
                raise Exception("Ax is given as a parameter. This function has no control over whether the figure "
                                + "will show.")
            plt.close(fig=fig)
        if save:
            if ax_given:
                raise Exception("Ax is given as a parameter. Please save outside of this function. \n" +
                                "Alternatively, try not giving ax as a parameter to this function")
            fig.savefig(f"{path}/{title}")
        if return_fig:
            # if ax_given:
            #     raise Exception("Ax is given as a parameter, and therefore fig
            return ax.figure

    def reference_with_xps(self, path="", scaling_factor=1):
        self.xps_df = xps.read_and_process(path, round=1) / scaling_factor

    def gaussian_df(self, t=-1, avg_duration=1, scaling_factor=1):
        if t == -1:
            t = self.df_raw.index.max()
        s = self.average_values(t, t + avg_duration) / scaling_factor

        # TODO
        # print(s)
        sigma = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)
        bes = [self.substances[name].orbitals[0].binding_energy for name in self.species_ordering]
        x_axis = np.arange(min(bes) - 5, max(bes) + 5, .1)
        envelope_vals = np.zeros(x_axis.size)

        res = {}

        for name in self.species_ordering:
            sol = s[name] if name in s else 0
            for o in self.substances[name].orbitals:
                be = o.binding_energy
                dist = sol * norm.pdf(x_axis, be, sigma)

                envelope_vals += dist
                res[name] = dist
        res["Envelope"] = envelope_vals
        res = pd.DataFrame.from_dict(res)
        res = res.set_index(x_axis)
        return res

    def bes(self):
        bes = {}
        for name in self.species_ordering:
            for o in self.substances[name].orbitals:
                bes[name] = o.binding_energy
        return bes

    @staticmethod
    def sample_for_xps(gaussian, xps_df):
        min_be = min(gaussian.index.min(), xps_df.index.min())
        max_be = max(gaussian.index.max(), xps_df.index.max())
        return xps.fill_zeros(gaussian, min_be, max_be), xps.fill_zeros(xps_df, min_be, max_be)

    def plot_gaussian(self, t=-1, avg_duration=1, path="", xps_path="", xps_scaling=1, save=False, return_fig=False,
                      fig_size="Default", dpi=100, scaling_factor=1, ax=None, envelope_name="CRN"):
        """
        Plot the Gaussian function from time t to time t + 1.
        """
        # if fig_size == "Default":
        #     plt.rcParams['figure.figsize'] = [27 / 2.54, 18 / 2.54]
        # elif fig_size == "Small":
        #     plt.rcParams['figure.figsize'] = [0.6 * 27 / 2.54, 0.6 * 18 / 2.54]
        # else:
        #     plt.rcParams['figure.figsize'] = fig_size
        gaussian = self.gaussian_df(t, avg_duration=avg_duration, scaling_factor=scaling_factor)
        min_be = gaussian.index.min()
        max_be = gaussian.index.max()
        bes = self.bes()

        curves = []
        if not ax:
            fig, ax = plt.subplots()
            if fig_size == "Default":
                fig.set_figheight(6)
                fig.set_figwidth(8)
            else:
                # For the plotting used in videos
                fig.set_dpi(dpi)
                fig.set_figheight(fig_size[1])
                fig.set_figwidth(fig_size[0])
        plt.style.use('seaborn-white')
        ax.tick_params(axis="x", direction="out", length=8, width=2)
        ax.tick_params(axis='both', which='minor', labelsize=12)

        gaussian_peaks = {}
        for n in sorted([n for n in gaussian.columns if n != "Envelope"],
                        key=lambda x: max(gaussian[x]), reverse=True):
            color = self.species_colors[n]
            peaks = ax.fill(gaussian.index, gaussian[n], label=n, color=color)
            gaussian_peaks[n] = peaks
            # Plot the Peak of the Gaussian
            ax.axvline(x=bes[n], ymin=0, ymax=1, color=color, linestyle='dashed', lw=1)

        legend_labels = []
        if xps_path:
            self.reference_with_xps(xps_path, xps_scaling)
            gaussian, xps_df = Results.sample_for_xps(gaussian, self.xps_df)
            min_be = min(gaussian.index.min(), xps_df.index.min())
            max_be = max(gaussian.index.max(), xps_df.index.max())

            curve = ax.plot(xps_df.index, xps_df["Envelope"], color='green', label="Experiment", alpha=0.5)
            curves.append(curve)
            legend_labels += ["Experiment"] + legend_labels

        curve = ax.plot(gaussian.index, gaussian["Envelope"], linewidth=2, color='black', label="CRN", alpha=0.8)
        curves.append(curve)
        legend_labels += [envelope_name] + [gaussian_peaks[n][0].get_label() for n in self.species_ordering]

        ax.legend([c[0] for c in curves] + [gaussian_peaks[n][0] for n in self.species_ordering], legend_labels,
                  fontsize=12)
        # if fig_size == "Default":
        #     ax.legend([c[0] for c in curves] + [gaussian_peaks[n][0] for n in names], legend_labels,
        #               fontsize=14)
        # elif fig_size == "Small":
        #     ax.legend([c[0] for c in curves] + [gaussian_peaks[n][0] for n in names], legend_labels,
        #               fontsize=12)

        ax.set_xlabel("Binding Energy", fontsize=12)
        ax.set_ylabel("Intensity", fontsize=12)

        ax.set_xlim(max_be, min_be)

        if save:
            ax.figure.savefig("{}/spectrum.png".format(path))
        if return_fig:
            # Usually, you don't need to return because plot would draw the graph in the current
            # Jupyter Notebook. Return only when you need the figure.
            return ax.figure

    def raw_string_gaussian(self, t=-1, avg_duration=1, y_upper_limit=None,  xps_scaling=1, fig_size="Default", dpi=100,
                            scaling_factor=1, ax=None):
        """
        A wrapper function for self.plot_gaussian intended for generating frames in Pygame videos.

        # TODO: explain the use of fig_size
        :return: raw string representation of the figure
        """
        backend = matplotlib.rcParams['backend']
        matplotlib.use("Agg")
        # TODO: determine if this solves the issue with inconsistent font.
        matplotlib.rcParams["font.family"] = "arial"
        fig = self.plot_gaussian(t=t, avg_duration=avg_duration, xps_scaling=xps_scaling, return_fig=True,
                                 fig_size=fig_size, dpi=dpi, scaling_factor=scaling_factor, ax=ax,
                                 envelope_name="total")
        fig.tight_layout()
        matplotlib.pyplot.ylim((0, y_upper_limit))
        canvas = agg.FigureCanvasAgg(fig)
        canvas.draw()
        renderer = canvas.get_renderer()
        matplotlib.use(backend)
        # TODO: close plt in case too many plots got opened.
        return renderer.tostring_rgb(),  canvas.get_width_height()

    def save(self, directory=None):
        """
        Save the data frame in a zipped file in directory.

        :param directory: a path in the file system; if None, infer the directory from the rules file.
        :return:
        """
        self.df_raw.to_csv("{}/Data_raw.gz".format(directory), compression='gzip')
        self.df.to_csv("{}/Data.gz".format(directory), compression='gzip')

    def save_converged_values(self, directory="", convergence_time=None, annotations=None):
        """
        Save converged values for fast access.

        :param directory:
        :param convergence_time: Time when convergence occurs; when not provided, provide an estimate.
        :param annotations: a dictionary of annotation name to annotation contents.
        :return:
        """
        if not convergence_time:
            convergence_time = self.estimate_convergence_time()
        s = self.compute_converged_values(convergence_time)
        species = list(s.index)
        values = [s[sps] for sps in species]
        results_succinct = {
            "Convergence Time": [convergence_time],
            "Species": species,
            "Final Values": values,
        }
        results_succinct.update(annotations)
        with open('{}/Results_Succinct.csv'.format(directory), 'w+') as f:
            writer = csv.writer(f)
            for k, v in results_succinct.items():
                writer.writerow([k] + v)

    def estimate_convergence_time(self):
        """
        :return: estimated  convergence time.
        """
        return self.df.index.max() - 20

    @property
    def last_time_stamp(self):
        return self.df_raw.index[-1]

    @property
    def max_concentration(self):
        return max([self.df_raw[s].max() for s in self.species_ordering])









