import numpy as np
from scipy.optimize import curve_fit

from lblcrn.common.curve import Curve
from lblcrn.common.util import unweave, weave
from lblcrn.xps_data_processing.line_shapes import line_shapes


def decompose(curve_df,
              function="glp",
              center_ranges=None,
              fwhm_ranges=None,
              mixing_param_means=None):
    """
    Decompose a curve into num_components independent Gaussian distributions, each of which
    has a mean that's in center_range and fwhm in fwhm_ranges.
    """
    f = produce_fitting_function(function)
    if center_ranges is not None and fwhm_ranges is not None:
        num_components = center_ranges.shape[0]

        # Translate full-width-half-maximum into standard deviation ranges.
        std_ranges = 2 * np.sqrt(2 * np.log(2)) * fwhm_ranges

        std_means = np.mean(std_ranges, axis=1)
        std_lower_bounds = std_ranges[:, 0]
        std_upper_bounds = std_ranges[:, 1]

        fwhm_means = np.mean(fwhm_ranges, axis=1)
        fwhm_lower_bounds = fwhm_ranges[:, 0]
        fwhm_upper_bounds = fwhm_ranges[:, 1]

        mean_means = np.mean(center_ranges, axis=1)
        mean_lower_bounds = center_ranges[:, 0]
        mean_upper_bounds = center_ranges[:, 1]

        # Give amplitude an upper bound of the maxima in assigned center ranges.
        amp_upper_bounds = np.array([curve_df[(curve_df.index >= mean_lower_bounds[i]) &
             (curve_df.index <= mean_upper_bounds[i])].iloc[:, 0].max() for i in range(num_components)])
        amp_means = 1/2 * amp_upper_bounds
        amp_lower_bounds = np.zeros(num_components)

        m_lower_bounds = np.zeros(num_components)
        m_upper_bounds = np.ones(num_components)
        if mixing_param_means is None:
            m_means = m_upper_bounds.copy()
        else:
            m_means = mixing_param_means

        initial_params = weave(amp_means.tolist(),  mean_means.tolist(), fwhm_means.tolist(), m_means.tolist())
        print("initial params", initial_params)

        lower_bound = weave(amp_lower_bounds.tolist(), mean_lower_bounds.tolist(), fwhm_lower_bounds.tolist(),
                            m_lower_bounds.tolist())
        upper_bound = weave(amp_upper_bounds.tolist(), mean_upper_bounds.tolist(), fwhm_upper_bounds.tolist(),
                            m_upper_bounds.tolist())

        print("bounds", [lower_bound, upper_bound])
        popt, pcov = curve_fit(f, np.flip(curve_df.index.to_numpy()), np.flip(curve_df.iloc[:, 0].to_numpy()),
                               p0=initial_params, bounds=[lower_bound, upper_bound])
    else:
        popt, pcov = curve_fit(f, np.flip(curve_df.index.to_numpy()), np.flip(curve_df.iloc[:, 0].to_numpy()))

    final_params = unweave(popt, line_shapes[function]["num_params"] - 1)

    print("result heights", final_params[0])
    print("result means", final_params[1])
    print("result widths", final_params[2])
    return Curve(f, curve_df.index.min(), curve_df.index.max(), len(curve_df.index),
                 curve_df, *final_params)


def produce_fitting_function(type="glp"):
    func = line_shapes[type]["function"]
    num_params = line_shapes[type]["num_params"] - 1

    def f(xs, *params):
        """
        Function of the assumed decomposition. This is the model; we optimize the params for this
        function towards certain goals.
        """
        params = unweave(params, num_params)
        ys = np.zeros_like(xs)
        for param_set in zip(*params):

            print("updated param_set", param_set)

            ys += func(xs, *param_set)
        return ys
    return f


def gaussian(xs, scale, loc, std):
    """
    A gaussian function for np arrays.
    """
    return scale * np.exp(- np.square(xs - loc) / (2 * std * std))
