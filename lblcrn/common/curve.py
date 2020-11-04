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
        self.f = function
        self.start, self.stop, self.n_samples = start, stop, n_samples
        self.params = params

        self.xs = np.flip(np.linspace(start, stop, n_samples))
        self.ys = np.zeros(n_samples)

        for p in zip(*params):
            self.ys += self.f(self.xs, *p)

        # for p in zip(*params):
        #     print(p)
        # print(self.ys)

        self.original_curve_df = original_curve_df

    def plot(self, target=None, data_fmt='b.', species_fmt='r--', verbose=True, data=False):
        if not target:
            target = plt.gca()
        # if data:  # Plot Data
        target.plot(self.xs, self.ys, data_fmt, label="Fitted Sum")

        for i, d in enumerate(zip(*self.params)):
            if verbose:
                label = f'Component {i}, params={[round(n, 2) for n in d]}'
            else:
                label = 'Component {}'.format(i)
            # Plot Species

            target.plot(self.xs, self.f(self.xs, *d), species_fmt, label=label)
        target.plot(self.original_curve_df.index, self.original_curve_df.iloc[:, 0], label="Original Curve")

        target.set_xlim(max(max(self.xs), self.original_curve_df.index.max()),
                        min(min(self.xs), self.original_curve_df.index.min()))
        target.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3), fancybox=True)

