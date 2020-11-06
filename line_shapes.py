"""
A collection of line shape functions.

The first 4 parameters for the functins are always:
xs: a numpy array of x values
h: height
e: center
f: width
"""

import numpy as np


def gaussian(x, h, e, f):
    return h * np.exp(-4 * np.log(2) * np.square(x - e) / (f ** 2))


def lorentzian(x, h, e, f):
    return h * 1 / (1 + 4*np.square((x - e) / f))


def la(x, h, e, f, alpha, beta, m):
    """
    Line shape signature: LA(alpha, beta, m)
    """
    y_lorentzian = np.vectorize(lambda a: lorentzian(a, h, e, f)**alpha if x < e else lorentzian(a, h, e, f)**beta)(x)
    y_gaussian = gaussian(x, h, e, m)
    return np.convolve(y_lorentzian, y_gaussian)


# Gaussian/Lorentzian Product Form
def glp(x, h, e, f, m):
    """
    Additional parameters:
    m: mixing parameter, ratio of Lorentzian.
    """
    return h * np.exp(-4 * np.log(2) * (1 - m) * np.square(x - e) / (f ** 2)) / (1 + 4*m*np.square(x - e)/(f**2))


# Gaussian/Lorentzian Product Form
def gls(x, h, e, f, m):
    return (1 - m) * gaussian(x, h, e, f) + m * lorentzian(x, h, e, f)


def tail_modifier(x, s, k, e, f):
    return s * np.exp(-k * (x - e) / f) if x <= e else 1


def glp_t(x, h, e, f, m, k):
    glp_result = glp(x, h, e, f, m)
    return glp_result + (1 - glp_result) * tail_modifier(x, 1, k, e, f)


def gls_t(x, h, e, f, m, k):
    gls_result = gls(x, h, e, f, m)
    return gls_result + (1 - gls_result) * tail_modifier(x, 1, k, e, f)


# Doniach-Sunjic profile
def doniach_sunjic(x, a, f, e):
    return np.math.gamma(a) * np.cos(np.pi * a / 2 + (1 - a) * np.arctan((x - e) / f)) / \
           (f**2 + np.square(x - e))**((1 - a) / 2)


# TODO: convolve with Gaussian
def doniach_sunjic_gaussian(x, a, f, e):
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
        "num_params": 6,
        "function": la
    },
    "ds": {
        "num"
    }
}

