"""Constants, mostly sympy-based, for the whole projec to depend on."""

from typing import Optional

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

def k_names(idx: Optional[int] = None):
    if idx is None: return 'k', 'κ'
    elif idx >= 0: return f'k_{idx}', f'κ_{idx}'

def k_map(idx: Optional[int] = None):
    k_forward, k_reverse = k_names()
    k_forward_dx, k_reverse_dx = k_names(idx)
    if idx is None:
        # With no index, it swaps forward and reverse
        return {
            k_forward: k_reverse_dx,
            k_reverse: k_forward_dx,
        }
    else:
        # With an index, it adds the
        return {
            k_forward: k_forward_dx,
            k_reverse: k_reverse_dx,
        }