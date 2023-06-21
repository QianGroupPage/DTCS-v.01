
from collections.abc import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython import display

from dtcs.twin.xps import XPSExperiment

import scipy
from datetime import datetime
from dtcs.optim.gp.core import kernel_l2_single_task
from gpcam.autonomous_experimenter import AutonomousExperimenterGP


def flatten_dictionary(dic, prefix=tuple()):
    # Base case
    if not isinstance(dic, dict):
        return {prefix: dic}

    out_dict = {}
    for key, value in dic.items():
        sub_dict = flatten_dictionary(value, prefix=prefix + (key,))
        out_dict.update(sub_dict)
    return out_dict


class CRNGibbsDataset:

    def __init__(self, sim, goal_concs=None, file_path=None, df=None):
        self.sim = sim
        self.goal_concs = {}
        self.df = None
        self.printer = self._default_printer
        self._fake_goal = None
        self._DEBUG_best_score = 1e10

        # Make the dataframe
        if file_path:
            # Load and save the data
            self.df = pd.read_csv(file_path, header=[0, 1], index_col=0)

            # Process the dataframe
            self._make_index()
        elif df is not None:
            self.df = df
        else:
            raise TypeError('Need to supply either file path or Dataframe.')

        # self.species = list(self.df['concs'].columns)
        # self.gibbs_cols = list(ds.gibbs.columns.droplevel(1).droplevel(1))
        self.set_goal(goal_concs if goal_concs is not None else {})

    def optimize_bh(
            dsg,
            fit: float = 1e-2,
            find: int = 1,
    ):
        # Time how long it takes to optimize
        start_time = datetime.now()

        # Make a function to handle local minima
        minima = []

        def local_minima_callback(gibbs, score, accepted):
            minima.append((gibbs, score, accepted))
            pass

        # Optimize
        try:
            scipy.optimize.basinhopping(
                func=lambda x: dsg.score(*x),
                x0=[0] * 7,
                minimizer_kwargs=dict(
                    bounds=[(-1, 1)] * 7, ),
                callback=local_minima_callback,
            )
        except KeyboardInterrupt:
            pass

    def optimize_gp(
            self,
            num_energies,  # TODO
    ):
        NUM_ENERGIES = num_energies

        def inst_func(to_score_list, *args, **kwargs):
            """Evaluate possible energies, using our scoring method.

            Args:
                to_score_list: A list of dicts. The 'position' key is the energies
                    that gpCAM wants us to test, and we set the 'value' key with
                    our scoring.

            Returns:
                A dictionary matching the structure of to_score_list, with the
                'value' keys filled in.
            """

            for to_score in to_score_list:
                gibbs = to_score['position']

                # Score the supplied energies
                score = self.score(*gibbs)
                to_score['value'] = score

            return to_score_list

        hpbounds = [[0.01, 1e5]] + [[0.001, 5000], ] * NUM_ENERGIES

        exp = AutonomousExperimenterGP(
            parameter_bounds=np.array([(-0.5, 0.5)] * NUM_ENERGIES),
            instrument_func=inst_func,
            hyperparameters=np.ones(NUM_ENERGIES + 1),
            hyperparameter_bounds=hpbounds,
            # acq_func='ucb',
            acq_func='covariance',
            x=[(0,) * NUM_ENERGIES],
            kernel_func=kernel_l2_single_task,
        )

        exp.train()
        exp.go(
            N=10000,
            acq_func_opt_setting=lambda _: 'global'  # lambda number: 'global' if number % 2 == 0 else 'global',
        )
        return exp

    @classmethod
    def from_dataset(cls, ds):
        """So I don't have to re-load the data each time I modify this class."""
        new_ds = cls(
            sim=ds.sim,
            goal_concs=ds.goal_concs,
            df=ds.df,
        )
        new_ds._fake_goal = ds._fake_goal
        return new_ds

    @classmethod
    def from_sim(cls, sim, num_energies, **kwargs):
        """Creates an empty dataset using a simulator function."""
        dic = sim([0] * num_energies)

        conc_dict = flatten_dictionary(dic, prefix=('concs',))
        mindex_len = len(next(iter(conc_dict)))
        gibbs_col_tupes = [('gibbs', f'G{ind}') + ('',) * (mindex_len - 2) for ind in range(1, num_energies + 1)]
        score_tupe = [('score',) + ('',) * (mindex_len - 1)]

        col_index = pd.Index(
            list(conc_dict.keys()) + gibbs_col_tupes + score_tupe
        )

        df = pd.DataFrame(
            columns=col_index,
        )
        return cls(
            sim=sim,
            df=df,
            **kwargs,
        )

    @property
    def scores(self):
        return self.df['score']

    @property
    def gibbs(self):
        return self.df['gibbs']

    @property
    def concs(self):
        return self.df['concs']

    def empty(self):
        self.df.drop(self.df.index, inplace=True)

    def fake_goal(self, goal_gibbs):
        """Makes fake goal data."""
        self.goal_concs = {}
        goal = self._simulate(goal_gibbs)
        goal['score'] = 0
        self._fake_goal = goal
        self.set_goal(self._fake_goal['concs'])
        return goal

    def set_goal(self, goal_concs):
        """Sets the goal concentrations and re-scores the dataframe."""
        self.goal_concs = goal_concs
        self.df['score'] = self._score(self.df['concs'])

    def score(self, *gibbs):
        """Check the score at the given Gibbs energies.

        If given new energies, will simulate and store the output."""
        # Make sure we get the right amount of energies
        # if len(gibbs) != len(self.gibbs):
        #     raise TypeError(f'Must supply exactly {len(self.gibbs)} energies.')

        # Check if we have that data already
        ridx = hash(gibbs)
        try:
            return self[ridx]['score'].iloc[0]
        except KeyError as err:
            self._sim_notarized(gibbs, ridx)
            return self[ridx]['score'].iloc[0]

    def plot(self, ridx):
        """Plot the given row."""
        # Simulate all the XPSs we need
        observables = []
        row = self.df.loc[ridx]['concs']

        # We're going to create title_and_concs, which is a list of
        #  names for each observable and the concentrations for it.
        #  It's more complicated if there's more than one condition.
        title_and_concs = []
        if row.index.nlevels > 1:
            # This groups by all but the last level in the series, which is the concentrations.
            for cond, series in row.groupby(level=tuple(range(row.index.nlevels))[:-1]):
                if not isinstance(cond, tuple):
                    cond = (cond,)
                concs = series[cond]
                title_and_concs.append((', '.join(cond), concs))
        else:
            # If there's only one set of conditions, it's easy
            title_and_concs = [('System', row)]

        for title, concs in title_and_concs:
            xps = XPSExperiment(
                species_manager=sm,
                # experimental=exp,
                sim_concs={specie: concs[specie.name] for specie in sm.symbols if specie.name in concs},
                autoresample=False,
            )

            observables.append((xps.resample(
                species=[o, oh, h2o, h2o_multi, o_h_h2o],
            ), title))

        # Make enough plots for all of them
        fig, axes = plt.subplots(
            len(observables), 1,
            figsize=(8, 2.5 * len(observables)),
            gridspec_kw=dict(
                hspace=0.5
            )
        )
        # If there's on plot, wrap it
        axes = axes if isinstance(axes, Iterable) else (axes,)

        for (xo, title), ax in zip(observables, axes):
            xo.plot(
                ax=ax,
                deconv_gaussians=False,
            )
            ax.set_title(title)

    @staticmethod
    def _default_printer(dsg, gibbs, ridx):
        display.clear_output(wait=True)
        gibbs_str = ', '.join(f'{gibb:.3f}' for gibb in gibbs)
        print(f'Simulating #{len(dsg.df)}')
        print(f'Row Index #{ridx}')
        print(f'Energies: {gibbs_str}')
        print(f'Current Best: {dsg._DEBUG_best_score}')

    def _sim_notarized(self, gibbs, ridx):
        """Simulates and records output."""
        self.printer(self, gibbs, ridx)

        row = self._simulate(gibbs)
        self.df.loc[ridx] = row

        score = float(row['score'])
        if score < self._DEBUG_best_score:
            self._DEBUG_best_score = score

    def _simulate(self, gibbs):
        """Simulates and scores the system with those energies."""
        concs_raw = flatten_dictionary(self.sim(gibbs))

        row = pd.Series(
            index=self.df.columns,
            dtype=np.float64,
        )

        for key, value in concs_raw.items():
            row[('concs',) + key] = value
        row['gibbs'] = gibbs
        row['score'] = self._score(row['concs'])
        return row

    def _score(self, concs):
        """Scores energies by via sum of log difference squared.

        If goal_concs is empty, will return 0."""
        score = 0
        for key, gconc in self.goal_concs.items():
            score += np.power(np.log(gconc) - np.log(concs[key]), 2)
        return score

    def _repr_html_(self):
        return self.df._repr_html_()

    def __getitem__(self, key):
        """Forwards to df.loc.__getitem__."""
        return self.df.loc[key]