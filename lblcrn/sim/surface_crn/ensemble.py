from lblcrn.sim.surface_crn.results import Results
import pandas as pd


class Ensemble:
    """
    Defines a solution to a system of differential equations, providing a variety of
    functions to manipulate and visualize the results
    """
    def __init__(self, results=None):
        self.results = results
        mean_df = pd.concat([r.df for r in self.results]).groupby(level=0).mean()
        self.mean_result = Results(results[0].manifest_file, results[0].rxns, mean_df)

    def plot_all(self):
        axes = Results.side_by_side_axes(len(self.results))

        x_max = max([r.last_time_stamp for r in self.results]) * 1.1
        y_max = max([r.max_concentration for r in self.results]) * 1.1

        for i, r in enumerate(self.results):
            ax = axes[i]
            r.plot_evolution(use_raw_data=True, ax=ax)
            ax.set_title(f"Result {i}")
            ax.set_xlim(0, x_max)
            ax.set_ylim(0, y_max)

    def plot_mean(self, ax=None, title=""):
        return self.mean_result.plot_evolution(use_raw_data=True, ax=ax, title=title)

    def __getitem__(self, arg):
        return self.results[arg]
