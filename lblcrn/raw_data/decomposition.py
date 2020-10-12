import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit
from lblcrn.common.gmm_fitting.gmm import Spectrum
from lblcrn.common.util import weave, unweave


def decompose(curve_df, center_ranges=None, fwhm_ranges=None):
    """
    Decompose a curve into num_components independent Gaussian distributions, each of which
    has a mean that's in center_range and fwhm in fwhm_ranges.
    """
    if center_ranges is not None and fwhm_ranges is not None:
        num_components = center_ranges.shape[0]

        # Translate full-width-half-maximum into standard deviation ranges.
        std_ranges = 2 * np.sqrt(2 * np.log(2)) * fwhm_ranges

        std_means = np.mean(std_ranges, axis=1)
        std_lower_bounds = std_ranges[:, 0]
        std_upper_bounds = std_ranges[:, 1]

        mean_means = np.mean(center_ranges, axis=1)
        mean_lower_bounds = center_ranges[:, 0]
        mean_upper_bounds = center_ranges[:, 1]

        # Give amplitude an upper bound of the maxima in assigned center ranges.
        amp_upper_bounds = np.array([curve_df[(curve_df.index >= mean_lower_bounds[i]) &
             (curve_df.index <= mean_upper_bounds[i])].iloc[:, 0].max() for i in range(num_components)])
        amp_means = 1/2 * amp_upper_bounds
        amp_lower_bounds = np.zeros(num_components)

        initial_params = weave(std_means.tolist(), mean_means.tolist(), amp_means.tolist())
        print("initial params", initial_params)

        lower_bound = weave(std_lower_bounds.tolist(), mean_lower_bounds.tolist(), amp_lower_bounds.tolist())
        upper_bound = weave(std_upper_bounds.tolist(), mean_upper_bounds.tolist(), amp_upper_bounds.tolist())
        popt, pcov = curve_fit(f, curve_df.index.to_numpy(), curve_df.iloc[:, 0].to_numpy(),
                            p0=initial_params)
    else:
        popt, pcov = curve_fit(f, curve_df.index.to_numpy(), curve_df.iloc[:, 0].to_numpy())

    stds, means, amps = unweave(popt, 3)

    print("result stds", stds)
    print("result means", means)
    print("result amps", amps)
    return Spectrum(stds, means, amps, curve_df.index.min(), 
                curve_df.index.max(), len(curve_df.index), original_curve_df=curve_df)


def f(xs, *params):
    """
    Function of the assumed decomposition. This is the model; we optimize the params for this 
    function towards certain goals.
    """
    stds, means, amps = unweave(params, 3)
    ys = np.zeros_like(xs)
    for std, mean, amp in zip(stds, means, amps):
        ys += gaussian(xs, amp, mean, std)
    return ys


def gaussian(xs, scale, loc, std):
    """
    A gaussian function for np arrays.
    """
    return scale * np.exp(- np.square(xs - loc)/ (2 * std * std))