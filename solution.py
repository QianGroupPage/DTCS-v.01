from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from typing import Callable, List, Dict, Tuple

class Solution:
    """
    Defines a solution to a system of differential equations, providing a variety of
    functions to manipulate and visualize the results
    """
    def __init__(self, t: List[float], y: List[List[float]], names: str):
        self.t = t
        self.states: Dict[str, List[float]] = {}
        self.wavelengths: Dict[str, float] = {}
        for i in range(len(y)):
            self.states[names[i]] = y[i]

    def sols(self) -> List[List[float]]:
        return self.states.values()

    def time_steps(self) -> List[float]:
        return self.t
    
    def basic_plot(self):
        """
        Draws a basic plot over the time domain
        """
        for sol in self.states.values():
            plt.plot(self.t, sol)

    def final_state(self) -> Dict[str, float]:
        """
        Returns the final state of the system of equations
        """
        final_state: Dict[str, float] = {}
        for name, sol in self.states.items():
            final_state[name] = sol[len(sol) - 1]
        return final_state

    def var_sols(self, *vars) -> Dict[str, List[float]]:
        """
        Returns a list of solutions for the specified variables passed in as arguments
        """
        return {name: sol for name, sol in self.states.items() if name in vars}

    def plot_gaussian(self):
        """
        Plots a gaussian distribution of the final species concentrations. FWHM is set at 0.75
        """
        sigma = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)
        
        for name, sol in self.final_state().items():
            wavelength = self.wavelengths[name]
            x_axis = np.arange(wavelength - 5, wavelength + 5, .001)
            plt.plot(x_axis, sol*norm.pdf(x_axis, wavelength, sigma))
        plt.show()
    
    def set_wavelengths(self, wavelengths: Dict[str, float]):
        """
        Assigns wavelength values to all species.
        """
        self.wavelengths = dict(wavelengths)

    def __repr__(self):
        shortened_sols = [sol[:10] for sol in self.states.values()]
        return 'Solution('+self.t[:10]+'..., '+shortened_sols[:10]+'...'


def solve_ode(ode: Callable[[float, List[float]], List[float]], names: List[str],
    time: float, init_vals: List[float], max_step: float = 0.1) -> Solution:
    """
    Solves a system of ordinary differential equations (described by a function) over a specified range of time.
    """
    sol = solve_ivp(ode, (0, time), init_vals, max_step=max_step)
    return Solution(sol.t, sol.y, names)
