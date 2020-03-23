from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
import numpy as np
from scipy.stats import norm
from typing import Callable, List, Dict, Tuple
import species
froem process_sympy_eqs import rxns_to_python_derivative_function, rxns_to_substances, rxns_to_initial_values

class Solution:
    """
    Defines a solution to a system of differential equations, providing a variety of
    functions to manipulate and visualize the results
    """
    def __init__(self, t: List[float], y: List[List[float]], rxns):
        # The time range of the reaction

        self.t = t
        # A map of symbols to concentrations over time
        self.states: Dict[str, List[float]] = {}
        # A map of symbols to substance information
        self.substances = dict(zip(rxns.get_symbols(), rxns.get_species()))
        #self.substances: Dict[str, species.Species] = {}
        for symbol, sub in self.substances.items():
        #    self.substances[substances[i].symbol] = substances[i]
            self.states[symbol] = y[rxns.symbol_index[symbol]]
        self.xps = None

        self.intensity = []

    def set_experimental(self, xps):
        """
        Takes an xps object contain experimental data, and scales the simulated solution
        to match.
        """
        self.xps = xps
        self.scaled_xps_intensity = []

    def sols(self) -> List[List[float]]:
        return list(self.states.values())

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

    def plot_gaussian(self, envelope: bool = False, overlay: bool = False):
        """
        Plots a gaussian distribution of the final species concentrations. FWHM is set at 0.75
        If specified, an envelope curve is also plotted
        """
        sigma = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)
        
        min_be = float('inf')
        max_be = float('-inf')
        
        # determine x axis bounds
        for substance in self.substances.values():
                min_be = min(substance.orbitals[0].binding_energy, min_be)
                max_be = max(substance.orbitals[0].binding_energy, max_be)

        x_axis = np.arange(min_be - 5, max_be + 5, .001)
        envelope_vals = np.zeros(x_axis.size)
        dists = []

        # plot a curve for each substance
        for name, sol in self.final_state().items():
            for o in self.substances[name].orbitals:
                be = o.binding_energy
                dist = sol * norm.pdf(x_axis, be, sigma)
                envelope_vals += dist
                dists.append(dist)

        for dist in sorted(dists, key=lambda x: max(x), reverse=True):
            plt.fill(x_axis, dist)

        if envelope:
            plt.plot(x_axis, envelope_vals, linewidth=4, color='black')
        
        if overlay:
            plt.plot(self.xps.binding_energy, self.xps.intensity, color='black')

        plt.gca().invert_xaxis()
        plt.show()

    def __repr__(self) -> str:
        shortened_sols = [sol[:10] for sol in self.states.values()]
        return 'Solution('+self.t[:10]+'..., '+shortened_sols[:10]+'...'


def solve_ode(ode: Callable[[float, List[float]], List[float]], rxns,
        time: float, rtol: float = 1e-3, atol: float = 1e-6, max_step: float = None) -> Solution:
    """
    Solves a system of ordinary differential equations (described by a function) over a specified range of time.

    Schedules should be constructed as follows: [[0, [1, 2]], [10, [3, 4]]], where the first number
    in the nested list is the time at which the changes should be made and the second list in the
    nested list is the changes to be made.
    """
    substances = rxns.get_species()

    # TODO: fix schedule!
    scheduleMap = {}
    for i, s in enumerate(substances):
        for t, state in s.schedule.schedule.items():
            if t not in scheduleMap:
                scheduleMap[t] = [0] * len(substances)
            scheduleMap[t][i] = state

    schedule = sorted(scheduleMap.items(), key=lambda x: x[0])
    #print(schedule)
    # ****

    concs = list(schedule[0][1])
    if len(schedule) == 1:
        if (max_step):
            sol = solve_ivp(ode, (0, time), concs, rtol=rtol, atol=atol, max_step=max_step)
        else:
            sol = solve_ivp(ode, (0, time), concs, rtol=rtol, atol=atol)

        y = sol.y
        t = sol.t
    else:
        t = np.empty(0)
        current_time = 0
        for s in schedule:
            for i, c in enumerate(s[1]):
                concs[i] += c

            sol = solve_ivp(ode, (current_time, s[0]), concs, rtol=rtol, atol=atol)
            if current_time == 0:
                t = sol.t
                y = sol.y
            else:
                t = np.append(t, sol.t)
                y = np.append(y, sol.y, axis=1)

            for i in range(len(concs)):
                concs[i] = y[i][len(y[i]) - 1]
            current_time = s[0]

    #print(y)
    return Solution(t, y, rxns)

def solve(rxns, time: float = 1, rtol: float = 1e-3, atol: float = 1e-6, max_step: float = None) -> Solution:
    return solve_ode(rxns_to_python_derivative_function(rxns), rxns, time,
            rtol=rtol, atol=atol, max_step=max_step)
