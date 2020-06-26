import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.stats import norm
import lblcrn.surface_crn.xps as xps
import os


class Results:
    """
    Defines a solution to a system of differential equations, providing a variety of
    functions to manipulate and visualize the results
    """
    TIME_COL_NAME = "Time (s) "  # column name used for the time entry

    def __init__(self, manifest_file, rxns, df=None):
        self.manifest_file = manifest_file
        self.df = df

        self.xps_df = None # The xps dataframe we have for reference

        self.species_ordering = None  # The ordering of species as they appear in our figures.
        self.species_colors = {}  # A dictionary from species names to their colors
        self.substances = dict(zip([repr(s) for s in rxns.get_symbols()], rxns.get_species())) if rxns else {}

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

    def average_values(self, start, end):
        """
        Average the result df from time start to time end.
        :return: a df with time column dropped
        """
        df = self.df
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
        for k, v in sub_species_dict.items():
            self.df[k] = pd.Series(0, self.df.index)
            for s in v:
                if s not in self.df.columns:
                    raise Exception(f"Subspecies {s} has not been recorded in the results data frame.")
                self.df[k] += self.df[s]

    # TODO: decrease the figure size in case zoom = False
    def plot_evolution(self, species_in_figure=None, start_time=0, end_time=-1, title="", save=False, path="",
                       zoom=False):
        """
        Plot the concentrations from start_time until time step end_time. -1 means till the end.

        Name is the specific name of this plot.
        """
        if end_time == -1:
            end_time = self.df.index.max()
        df = self.df[(self.df.index <= end_time) & (self.df.index >= start_time)]

        if species_in_figure is None:
            species_in_figure = self.species_ordering

        fig = plt.figure(figsize=(30, 60))
        fig.subplots_adjust(top=0.95)
        fig.suptitle(f"{title}", fontsize=48)

        if zoom:
            irange = range(len(species_in_figure) - 1)
        else:
            irange = [0]
        for i in irange:
            plt.subplot(len(species_in_figure), 1, i + 1)
            for j in range(i, len(species_in_figure)):
                species = species_in_figure[j]
                plt.tick_params(axis='both', which='both', labelsize=36)
                plt.plot(df[species], color=self.species_colors[species], label=species, linewidth=5)
                plt.legend(fontsize=36, numpoints=30)

            plt.xlabel("Time (s)", fontsize=36)
            plt.ylabel("Molecule Count (#)", fontsize=36)
        if save:
            fig.savefig(f"{path}/{title}")
        plt.close()
        return fig

    def reference_with_xps(self, path="", scaling_factor=1):
        self.xps_df = xps.read_and_process(path, round=1) / scaling_factor

    def gaussian_df(self, t=-1, scaling_factor=1):
        if t == -1:
            t = self.df.index.max()
        s = self.average_values(t, t + 1) / scaling_factor

        sigma = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)
        bes = [self.substances[name].orbitals[0].binding_energy for name in self.species_ordering]
        x_axis = np.arange(min(bes) - 5, max(bes) + 5, .1)
        envelope_vals = np.zeros(x_axis.size)

        res = {}

        for name in self.species_ordering:
            sol = s[name]
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

    def plot_gaussian(self, t=-1, path="", xps_path="", xps_scaling=1, save=False, show_fig=True, fig_size="Default",
                      scaling_factor=1,
                      ax=None):
        """
        Plot the Gaussian function from time t to time t + 1.
        """
        # if fig_size == "Default":
        #     plt.rcParams['figure.figsize'] = [27 / 2.54, 18 / 2.54]
        # elif fig_size == "Small":
        #     plt.rcParams['figure.figsize'] = [0.6 * 27 / 2.54, 0.6 * 18 / 2.54]
        # else:
        #     plt.rcParams['figure.figsize'] = fig_size
        gaussian = self.gaussian_df(t, scaling_factor)
        min_be = gaussian.index.min()
        max_be = gaussian.index.max()
        bes = self.bes()

        curves = []
        if not ax:
            fig, ax = plt.subplots()
        plt.style.use('seaborn-white')
        # ax.tick_params(axis="x", direction="out", length=8, width=2)
        # plt.tick_params(axis='both', which='minor', labelsize=36)

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

        curve = ax.plot(gaussian.index, gaussian["Envelope"], linewidth=4, color='black', label="CRN", alpha=0.8)
        curves.append(curve)
        legend_labels += ["CRN"] + [gaussian_peaks[n][0].get_label() for n in self.species_ordering]

        ax.legend([c[0] for c in curves] + [gaussian_peaks[n][0] for n in self.species_ordering], legend_labels)
        # if fig_size == "Default":
        #     ax.legend([c[0] for c in curves] + [gaussian_peaks[n][0] for n in names], legend_labels,
        #               fontsize=14)
        # elif fig_size == "Small":
        #     ax.legend([c[0] for c in curves] + [gaussian_peaks[n][0] for n in names], legend_labels,
        #               fontsize=12)

        ax.set_xlim(max_be, min_be)
        if show_fig:
            plt.show()
        if save:
            ax.figure.savefig("{}/spectrum.png".format(path))
        return ax.figure

    def save(self, directory=None):
        """
        Save the data frame in a zipped file in directory.

        :param directory: a path in the file system; if None, infer the directory from the rules file.
        :return:
        """
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







