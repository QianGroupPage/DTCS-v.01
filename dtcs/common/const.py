"""Constants, mostly sympy-based, for the whole projec to depend on."""

from sympy.physics import units
import sympy as sym

T = TIME = sym.Symbol('\mathbf{t}')
P = PRESSURE = sym.Symbol('\mathbf{P}')
K = TEMPERATURE = sym.Symbol('\mathbf{T}')
DG = GIBBS_ENERGY = sym.Symbol('\mathbf{\Delta G}')

RESERVED_SYMBOLS = [T, P, K, DG]

PRETTY_SUBS = {
    GIBBS_ENERGY: sym.Symbol('\Delta G'),
    PRESSURE: sym.Symbol('P'),
    TEMPERATURE: sym.Symbol('T'),
    TIME: sym.Symbol('t'),
    units.boltzmann: sym.Symbol('k_B')
}