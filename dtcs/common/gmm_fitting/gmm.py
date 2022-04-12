import numpy as np
from dtcs.common.util import weave, unweave
from scipy.stats import norm, uniform
from scipy.optimize import curve_fit


class XPSModel:
    def __init__(self, std=uniform(loc=.7, scale=0),
                 mean=uniform(loc=5, scale=10),
                 amp=uniform(loc=1, scale=9),
                 n=5, start=2.5, stop=17.5, n_samples=500):
        self.n, self.std, self.mean, self.amp = n, std, mean, amp
        self.start, self.stop, self.n_samples = start, stop, n_samples

    def _SM(self, stds, means, amps):
        return Spectrum(stds, means, amps, self.start, self.stop, self.n_samples)

    def sample(self):
        stds = self.std.rvs(self.n)
        means = self.mean.rvs(self.n)
        amps = self.amp.rvs(self.n)
        return self._SM(stds, means, amps)

    def fit(self, target_model, guess_n=None):
        if not guess_n:
            guess_n = self.n

        # Initial Params
#         stds = [self.std.mean()] * guess_n
        stds = [.7] * guess_n
        means = np.linspace(self.mean.ppf(.05), self.mean.ppf(.95), guess_n)
        amps = [self.amp.mean()] * guess_n
        initial_params = weave(stds, means, amps)
        # Lower Bound
#         stds = [self.std.ppf(0)] * guess_n
        stds = [.699999] * guess_n
        means = [self.mean.ppf(0)] * guess_n
        amps = [0] * guess_n
        lower_bound = weave(stds, means, amps)
        # Upper Bound
#         stds = [self.std.ppf(1)] * guess_n
        stds = [.7] * guess_n
        means = [self.mean.ppf(1)] * guess_n
        amps = [self.amp.ppf(1)] * guess_n
        upper_bound = weave(stds, means, amps)

        popt, pcov = curve_fit(Spectrum.f, target_model.xs, target_model.ys,
                               p0=initial_params, bounds=(lower_bound, upper_bound))
        return self._SM(*unweave(popt, 3))


class Spectrum:
    id = 1

    def __init__(self, stds, means, amps,
                 start, stop, n_samples, eps=.1):
        self.stds, self.means, self.amps = stds, means, amps
        self.start, self.stop, self.n_samples = start, stop, n_samples
        self.params = weave(stds, means, amps)

        # Sample
        self.xs = np.linspace(start, stop, n_samples)
        self.ys = self.pdf(self.xs) + norm.rvs(size=len(self.xs), scale=eps)

        self.id = Spectrum.id
        Spectrum.id += 1
        self.name = "Spectrum {}".format(self.id)

    @classmethod
    def f(cls, xs, *params):
        stds, means, amps = unweave(params, 3)
        ys = np.zeros_like(xs)
        for std, mean, amp in zip(stds, means, amps):
            ys += norm.pdf(xs, loc=mean, scale=std) * amp
        return ys

    def pdf(self, xs):
        ys = np.zeros_like(xs)
        for std, mean, amp in zip(self.stds, self.means, self.amps):
            ys += norm.pdf(xs, loc=mean, scale=std) * amp
        return ys

    def plot(self, target, data_fmt='b.', species_fmt='r--', verbose=True, data=False):
        if data:  # Plot Data
            target.plot(self.xs, self.ys, data_fmt, label='{} Data'.format(self.name))

        for i, d in enumerate(zip(self.stds, self.means, self.amps)):
            std, mean, amp = d
            if verbose:
                label = '{} Species {}: std={:.3} mean={:.3} amp={:.3}'.format(self.name, i, std, mean, amp)
            else:
                label = '{} Species {}'.format(self.name, i)
            # Plot Species
            target.plot(self.xs, norm.pdf(self.xs, loc=mean, scale=std) * amp, species_fmt, label=label)
