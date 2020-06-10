from scipy.integrate import solve_ivp, trapz
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, mean_absolute_error
from math import sqrt

from matplotlib.widgets import CheckButtons
import numpy as np
from scipy.stats import norm
from typing import Callable, List, Dict, Tuple

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

        self.envelope = []
        self.binding_energies = []
        self.resampled_intensity = []
        self.distributions = []
        self.ignore = []
        self.xps = None

    def set_experimental(self, xps):
        """Takes an xps object contain experimental data, and scales the simulated solution.
        """
        self.xps = xps

    def sols(self) -> List[List[float]]:
        return list(self.states.values())

    def time_steps(self) -> List[float]:
        return self.t
    
    def basic_plot(self):
        """Draws a basic plot over the time domain
        """
        for name, sol in self.states.items():
            plt.plot(self.t, sol, label=name)
        plt.legend()

    def final_state(self) -> Dict[str, float]:
        """Returns the final state of the system of equations
        """
        final_state: Dict[str, float] = {}
        for name, sol in self.states.items():
            final_state[name] = sol[len(sol) - 1]
        return final_state

    def var_sols(self, *vars) -> Dict[str, List[float]]:
        """Returns a list of solutions for the specified variables passed in as arguments
        """
        return {name: sol for name, sol in self.states.items() if name in vars}
    
    def set_ignore(self, ignore):
        """
        Set species to be ignored
        """
        self.ignore = list(ignore)

    def scale(self, to_scale, exp):
        """Scale experimental data intensity to match that of the experimental data.

        The largest value of the simulated data is scaled to match the largest value of the
        experimental data, and this scaling factor is then applied across all simulated data.
        """
        max_to_scale = max(to_scale[1])
        max_exp = max(exp[1])
        new_envelope = []
        new_dists = []
        scaling = max_exp / max_to_scale

        for v in to_scale[1]:
            new_envelope.append(v * scaling)

        for d in to_scale[2]:
            new_dists.append(d * scaling)

        return new_envelope, new_dists

    def resample(self):
        """Resample the simulated data, reducing its size to match that of the experimental data.
        """
        # rei = []
        # i = 0
        # bes = self.binding_energies
        # for intensity, be in zip(list(reversed(self.xps.intensity)), list(reversed(self.xps.binding_energy))):
            # while i < len(bes) and be > bes[i]:
                # rei.append(intensity)
                # i += 1
            # if i >= len(bes):
                # break
        # self.resampled_intensity = rei
        e = self.envelope
        bes = self.binding_energies
        xbes = self.xps.binding_energy
        xi = self.xps.intensity

        r_e, r_bes = [], []
        i = 0
        for b in list(reversed(xbes)):
            while  i < len(e) and bes[i] < b:
                i += 1
            r_e.append(e[i])
            r_bes.append(b)
        self.envelope = r_e
        self.resampled_binding_energies = np.array(r_bes)
        self.resampled_intensity = list(reversed(xi))

    def process(self, gas_range=()): 
        """Resample, scale, and calculate envelopes and other characteristics of the data.

        First the binding energy bounds are found, the envelope curve and individual species
        gaussians are then computed. Finally, if an experimental data object has been set, the
        simulated data is resampled and sacled.
        """
        sigma = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)
        
        min_be = float('inf')
        max_be = float('-inf')
        
        # determine x axis bounds
        for name, substance in self.substances.items():
            if name not in self.ignore:
                min_be = min(substance.orbitals[0].binding_energy, min_be)
                max_be = max(substance.orbitals[0].binding_energy, max_be)

        self.binding_energies = np.arange(min_be - 5, max_be + 5, .001)
        self.envelope = np.zeros(self.binding_energies.size)
        self.distributions = []
        self.names = []

        for name, sol in self.final_state().items():
            if name not in self.ignore:
                for o in self.substances[name].orbitals:
                    be = o.binding_energy
                    dist = sol * norm.pdf(self.binding_energies, be, sigma)
                    self.envelope += dist
                    self.distributions.append(dist)
                    self.names.append(name)

        self.resampled_binding_energies = self.binding_energies

        if self.xps:
            self.resample()
            self.envelope, self.distributions = self.scale((self.resampled_binding_energies, self.envelope,
                self.distributions), (self.xps.binding_energy, self.xps.intensity))

            # If a valid gas range is specified, create a fake peak
            if len(gas_range) == 2:
                start = gas_range[0]
                end = gas_range[1]
                intensities = list(reversed(self.xps.intensity))
                bes = list(reversed(self.xps.binding_energy))

                i = 0
                while i < len(bes) and bes[i] < start:
                    i += 1

                if i < len(bes):
                    max_intensity = intensities[i]
                    be = i
                    while i < len(bes) and i < bes[i] < end:
                        if intensities[i] > max_intensity:
                            max_intensity = intensities[i]
                            be = bes[i]
                        i += 1
                    if i < len(bes):
                        dist = max_intensity * norm.pdf(self.binding_energies, be, sigma)
                        resampled_dist = max_intensity * norm.pdf(self.resampled_binding_energies, be, sigma)
                        self.envelope += resampled_dist
                        self.distributions.append(dist)
                        self.names.append('gas phase')

    def plot_gaussian(self, envelope: bool = False, overlay: bool = False, resample_envelope: bool =
            False, ax=None, title=''):
        """
        Plots a gaussian distribution of the final species concentrations. FWHM is set at 0.75
        If specified, an envelope curve is also plotted
        """
        colors = ['red', 'green', 'orange', 'blue', 'purple', 'pink', 'yellow', 'gray', 'cyan']
        if not ax:
            for i, dist in sorted(enumerate(self.distributions), key=lambda x: max(x[1]), reverse=True):
                plt.fill(self.binding_energies, dist, label=self.names[i], color=colors[i])
            plt.legend()

            if overlay:
                if resample_envelope:
                    plt.plot(self.resampled_binding_energies, self.resampled_intensity, color='green')
                else:
                    plt.plot(self.xps.binding_energy, self.xps.intensity, color='green')

            if envelope:
                plt.plot(self.resampled_binding_energies, self.envelope, linewidth=4, color='black')
            
            plt.gca().invert_xaxis()
            plt.show()
        else:
            for i, dist in sorted(enumerate(self.distributions), key=lambda x: max(x[1]), reverse=True):
                ax.fill(self.binding_energies, dist, label=self.names[i], color=colors[i])
            ax.legend()

            if overlay:
                if resample_envelope:
                    ax.plot(self.resampled_binding_energies, self.resampled_intensity, color='green')
                else:
                    ax.plot(self.xps.binding_energy, self.xps.intensity, color='green')

            if envelope:
                ax.plot(self.resampled_binding_energies, self.envelope, linewidth=4, color='black')
            
            ax.set_xlim(max(self.resampled_binding_energies), min(self.resampled_binding_energies))
            ax.title.set_text(title)

    def rmse(self):
        return sqrt(mean_squared_error(self.resampled_intensity, self.envelope))

    def mae(self):
        return mean_absolute_error(self.resampled_intensity, self.envelope)
    
    def integral_diff(self):
        return abs(trapz(self.resampled_intensity, self.resampled_binding_energies) -
                trapz(self.envelope, self.resampled_binding_energies))

    def __repr__(self) -> str:
        shortened_sols = [sol[:10] for sol in self.states.values()]
        return 'Solution('+str(self.t[:10])+'..., '+str(shortened_sols[:10])+'...'