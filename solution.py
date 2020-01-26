from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
import numpy as np
from scipy.stats import norm
from typing import Callable, List, Dict, Tuple
import substance

class Solution:
    """
    Defines a solution to a system of differential equations, providing a variety of
    functions to manipulate and visualize the results
    """
    def __init__(self, t: List[float], y: List[List[float]], substances: List[substance.Substance]):
        # The time range of the reaction
        self.t = t
        # A map of symbols to concentrations over time
        self.states: Dict[str, List[float]] = {}
        # A map of symbols to substance information
        self.substances: Dict[str, substance.Substance] = {}
        for i in range(len(y)):
            self.substances[substances[i].symbol] = substances[i]
            self.states[substances[i].symbol] = y[i]

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

    def plot_gaussian(self, show_envelope: bool = False):
        """
        Plots a gaussian distribution of the final species concentrations. FWHM is set at 0.75
        If specified, an envelope curve is also plotted
        """
        sigma = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)
        
        min_binding_energy = float('inf')
        max_binding_energy = float('-inf')
        
        # determine x axis bounds
        for substance in self.substances.values():
            if substance.binding_energy < min_binding_energy:
                min_binding_energy = substance.binding_energy
            if substance.binding_energy > max_binding_energy:
                max_binding_energy = substance.binding_energy
        x_axis = np.arange(min_binding_energy- 5, max_binding_energy + 5, .001)

        envelope = np.zeros(x_axis.size)
        # plot a curve for each substance
        for name, sol in self.final_state().items():
            binding_energy = self.substances[name].binding_energy
            distribution = sol * norm.pdf(x_axis, binding_energy, sigma)
            envelope += distribution
            plt.plot(x_axis, distribution)

        if show_envelope:
            plt.plot(x_axis, envelope)

        plt.show()

    def __repr__(self):
        shortened_sols = [sol[:10] for sol in self.states.values()]
        return 'Solution('+self.t[:10]+'..., '+shortened_sols[:10]+'...'


def solve_ode(ode: Callable[[float, List[float]], List[float]], substances: List[substance.Substance],
    time: float, init_vals: List[float], max_step: float = 0.1) -> Solution:
    """
    Solves a system of ordinary differential equations (described by a function) over a specified range of time.
    """
    sol = solve_ivp(ode, (0, time), init_vals, max_step=max_step)
    return Solution(sol.t, sol.y, substances)
