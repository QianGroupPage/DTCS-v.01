"""
A collection of line shape functions.

The first 4 parameters for the functions are always:
xs: a numpy array of x values
height: height
center: center
fwhm: width
"""

import numpy as np

# --- Basic Line Shapes ------------------------------------------------\


def gaussian(x, height, center, fwhm):
    return height * np.exp(-4 * np.log(2) * np.square(x - center) / (fwhm ** 2))


def lorentzian(x, height, center, fwhm):
    return height * 1 / (1 + 4*np.square((x - center) / fwhm))


def tail_modifier(x, s, k, center, fwhm):
    return s * np.exp(-k * (x - center) / fwhm) if x <= center else 1


# --- Mixed Shapes ------------------------------------------------\


def glp(x, height, center, fwhm, m):
    """
    Gaussian/Lorentzian Product Form

    Line shape signature: GL(m).
    Additional parameter:
    m: mixing parameter, ratio of Lorentzian.
    """
    return height * np.exp(-4 * np.log(2) * (1 - m) * np.square(x - center) /
                           (fwhm ** 2)) / (1 + 4*m*np.square(x - center)/(fwhm**2))


def la(x, height, center, fwhm, alpha, beta, m):
    """
    Lorentzian Assymetric Form

    Line shape signature: LA(alpha, beta, m).
    Additional parameters:
    alpha: power applied to the Lorentzian in the region left of center; skews the left half;
    beta: power applied to the Lorentzian in the region right of center; skews the right half;
    m: fwhm of the Gaussian used to convolve with the skewed Lorentzian, a percentage value in the range of (0, 500).
    """
    m /= 100
    y_lorentzian = np.vectorize(lambda a: np.power(lorentzian(a, 1, center, fwhm), alpha) if a < center else
                   np.power(lorentzian(a, 1, center, fwhm), beta))(x)
    y_gaussian = gaussian(x, 1, center, m)
    return height * np.convolve(y_lorentzian, y_gaussian, mode="same")


# --- Less-used Shapes ------------------------------------------------\


def gls(x, height, center, fwhm, m):
    """
    Gaussian/Lorentzian Sum Form

    Line shape signature: GLS(m).
    Additional parameter:
    m: mixing parameter, ratio of Lorentzian.
    """
    return (1 - m) * gaussian(x, height, center, fwhm) + m * lorentzian(x, height, center, fwhm)


def glp_t(x, h, e, f, m, k):
    glp_result = glp(x, h, e, f, m)
    return glp_result + (1 - glp_result) * tail_modifier(x, 1, k, e, f)


def gls_t(x, h, e, f, m, k):
    gls_result = gls(x, h, e, f, m)
    return gls_result + (1 - gls_result) * tail_modifier(x, 1, k, e, f)


def doniach_sunjic(x, a, f, e):
    """
    Doniach-Sunjic profile
    """
    return np.math.gamma(a) * np.cos(np.pi * a / 2 + (1 - a) * np.arctan((x - e) / f)) / \
           (f**2 + np.square(x - e))**((1 - a) / 2)


def doniach_sunjic_gaussian(x, a, f, e):
    # TODO: convolve with Gaussian
    pass


line_shapes = {
    "glp": {
        "num_params": 5,
        "function": glp
    },
    "gls": {
        "num_params": 5,
        "function": gls
    },
    "glp_t": {
        "num_params": 6,
        "function": glp_t
    },
    "gls_t": {
        "num_params": 6,
        "function": gls_t
    },
    "la": {
        "num_params": 7,
        "function": la
    },
    "ds": {
        "num_params": 7,
        "function": doniach_sunjic
    }
}
