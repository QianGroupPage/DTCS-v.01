"""
Wrapper and visualizer for a curve based on a function and its list of parameter values.

Created by Ye Wang, based on work by Domas Buracas.
"""

import matplotlib.pyplot as plt
import numpy as np


class Curve:
    def __init__(self, function, start, stop, n_samples,
                 original_curve_df,
                 *params):
        """

        :param function:
        :param start:
        :param stop:
        :param n_samples:
        :param original_curve_df:
        :param params:
        """
        self.f = function
        self.start, self.stop, self.n_samples = start, stop, n_samples
        self.params = params

        self.xs = np.flip(np.linspace(start, stop, n_samples))
        self.ys = np.zeros(n_samples)

        for p in zip(*params):
            print("parameter  set", p)
            self.ys += self.f(self.xs, *p)

        self.original_curve_df = original_curve_df

    def plot(self, ax=None, data_fmt='b.', species_fmt='r--', verbose=True, data=False):
        if not ax:
            ax= plt.gca()
        # if data:  # Plot Data
        ax.plot(self.xs, self.ys, data_fmt, label="Fitted Sum")

        for i, d in enumerate(zip(*self.params)):
            if verbose:
                label = f'Component {i}, params={[round(n, 2) for n in d]}'
            else:
                label = 'Component {}'.format(i)
            # Plot Species
            ax.plot(self.xs, self.f(self.xs, *d), species_fmt, label=label)
        ax.plot(self.original_curve_df.index, self.original_curve_df.iloc[:, 0], label="Original Curve")

        ax.set_xlim(max(max(self.xs), self.original_curve_df.index.max()),
                        min(min(self.xs), self.original_curve_df.index.min()))
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3), fancybox=True)
