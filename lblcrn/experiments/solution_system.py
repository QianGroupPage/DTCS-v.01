from typing import List, Dict
import matplotlib.pyplot as plt
import sympy as sym
from lblcrn.crn_sym import reaction
import lblcrn.experiments.xps_io as xps_io
import lblcrn.experiments.xps as xps
import pandas as pd

class SolutionSystem:
    """Process and visualize a system of experiment

    This system performs automatic global scaling as well as other processing.

    Instance Variables:
        ignore: A list of sympy symbols to ignore
        systems: A list of diction
    """
    def __init__(self, experiments: List, experimental_files: List[str]):
        """Create a new solution system

        Given a list of lists of XPSExperiments (for all variations that are simulated) and a
        corresponding list of experimental files, experimental data is added to the XPSExperiment
        objects and the data is then auto-scaled.
        """
        self.systems = experiments

        # Iterate through XPS classes, adding experimental data.
        for i, f in enumerate(experimental_files):
            xps_exp = xps_io.read_new_data(f)[0]
            series = pd.Series(data=xps_exp.intensities, index=xps_exp.binding_energies)
            for s in self.systems[i]:
                s.experimental = series

        scaling_factor, max_index = self.max_scaling_factor()
        self.scale(scaling_factor)
        self._ignore = []

        print('scaling factor:', scaling_factor, '\tmax index:', max_index)

    def _process(self):
        """Rescale solution data.
        """
        scaling_factor, max_index = self.max_scaling_factor()
        self.scale(scaling_factor)

    @property
    def ignore(self) -> List[sym.Symbol]:
        return self._ignore
    
    @ignore.setter
    def ignore(self, ignore: List[sym.Symbol]) -> None:
        self._ignore = ignore
        for sys in self.systems:
            for s in sys:
                s.ignore = self._ignore
        self._process()

    def max_scaling_factor(self) -> int:
        """Find the maximum peak and return the scaling factor for that data
        and the system index."""
        self.systems[0][0].scale_factor = 1
        max_env = max(self.systems[0][0].envelope)
        max_index = -1
        max_exp = max(self.systems[0][0].experimental)

        # TODO: I changed this as little as possible but it seems to put
        # TODO: envelope to be larger than peak.
        # Find max peak of systems
        for i in range(len(self.systems)):
            for j in range(len(self.systems[i])):
                self.systems[i][j].scale_factor = 1
                potential = max(self.systems[i][j].experimental)
                if potential > max_exp:
                    max_exp = potential
                    max_index = i
                    max_env = max(self.systems[i][j].envelope)

        return max_exp / max_env, max_index

    def scale(self, scaling_factor):
        """Scale experimental data intensity to match that of the experimental data.
        """
        for sys in self.systems:
            for s in sys:
                s.scale_factor = scaling_factor
    
    def __getitem__(self, index):
        return self.systems[index]
    
    def plot(self, index, rows=0, cols=0):
        fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=(60, 60))
        plt.figure(fig.number)

        for j, ax in enumerate(fig.axes):
            plt.sca(ax)
            self.systems[index][j].plot(show=False)
            # TODO(Rithvik) I might have broken this line, sorry.
            plt.title(f'Eq: {str(int(j / rows))} Const: {str(j % cols)}')
        plt.show()

class XPSInitializationData:
    """Stores the reaction constants and raw experimental file location.

    This class is simply syntactic sugar for a map.

    Instance Variables:
        name: The name of this specific experiment.
        constants: A list reaction constants for this experiment, with an order corresponding to
            that of the reaction system to be simulated.
        experimental_file: The location of the experimental XPS file to be used for comparison
    """
    def __init__(self, title: str, constants: List[float], experimental_file: str):
        self.title = title
        self.constants = constants
        self.experimental_file = experimental_file


class XPSSystemRunner:
    """Orchestrate an end-to-end XPS analysis pipeline.
    
    Given a reaction system and initialization data objects, runs the specified simulations, returning
    a solution system object that can be used to visualize/analyze results.

    Instance Variables:
        rsys_generator: A function that takes in a list of reaction constants and returns a reaction
        system to be used.
        time: The amount of time for which the simulation is run.
        initialization_data: A list of XPSInitializationData objects that specify the various
        simulations to be run.
        multipliers: A list of reaction constant multipliers to be used.
        solutions: A list of lists of simulated data for an XPS experiment.
        complete: A list of boolean values representing whether all systems have been simulated.
    """

    def __init__(self,
                 rsys_generator,
                 time: float,
                 initializer_data: List[XPSInitializationData],
                 multipliers: List[float]):
        """Create a new XPS system runner.
        """
        self.rsys_generator = rsys_generator
        self.time = time
        self.initializer_data = initializer_data
        self.multipliers = multipliers
        self.solutions = [None for _ in range(len(self.initializer_data))]
        self.complete = [False for _ in range(len(self.initializer_data))]

    def simulate(self, index: int):
        # Test code that will be moved into a module
        data = self.initializer_data[index]
        sols = []

        for i in range(len(data.constants)):
            for j in range(len(self.multipliers)):
                scaled = list(data.constants)
                scaled[i] *= self.multipliers[j]

                rsys = self.rsys_generator(scaled)
                s = xps.simulate_xps(rsys, time=self.time, title=data.title + " Eq: "
                    + str(i) + "Constant: " + str(j))

                sols.append(s)
                print('Solved for ('+str(i)+', '+str(j)+')')
                print(scaled)
                print('\n')
        self.solutions[index] = sols
        self.complete[index] = True
                
    def system(self) -> SolutionSystem:
        """Returns a solution system if all experiments have been simulated.
        """
        if not all(self.complete):
            print("All experiments have not yet been simulated.")
            return
        return SolutionSystem(self.solutions, [x.experimental_file for x in self.initializer_data])
