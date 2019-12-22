from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from typing import Callable, List

class Solution:
    """
    Defines a solution to a system of differential equations, providing a variety of
    functions to manipulate and visualize the results
    """
    def __init__(self, t: List[float], y: List[List[float]]):
        self.t = t
        self.y = y

    def sols(self) -> List[List[float]]:
        return self.y

    def time_steps(self) -> List[float]:
        return self.t
    
    def basic_plot(self):
        """
        Draws a basic plot over the time domain
        """
        for sol in self.y:
            plt.plot(self.t, sol)

    def final_state(self) -> List[float]:
        """
        Returns the final state of the system of equations
        """
        final_state: List[float] = []
        for sol in self.y:
            final_state.append(sol[len(sol) - 1])
        return final_state

    def var_sols(self, *vars) -> List[List[float]]:
        """
        Returns a list of solutions for the specified variables passed in as arguments
        """
        return [sol for i, sol in enumerate(self.y) if i in vars]

    def plot_gaussian(self):
        final_state = self.final_state()
        sigma = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)
        
        #TODO(rithvik): take in wavelengths
        x_axis = np.arange(-5, len(final_state)+5, 0.001)
        for i, s in enumerate(final_state):
            plt.plot(x_axis, s*norm.pdf(x_axis, i, sigma))
        plt.show()

    def __repr__(self):
        shortened_y = [sol[:10] for sol in self.y]
        return 'Solution('+self.t[:10]+'..., '+shortened_y[:10]+'...'


def solve_ode(odes: List[Callable[[float, List[float]], List[float]]],
    time: float, init_vals: List[float], max_step: float = 0.1) -> Solution:
    """
    Solves a system of ordinary differential equations (described by functions) over a specified range of time.
    """
    sol = solve_ivp(odes, (0, time), init_vals, max_step=max_step)
    return Solution(sol.t, sol.y)
