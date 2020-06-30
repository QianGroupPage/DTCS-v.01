from typing import List
import matplotlib.pyplot as plt
import sympy as sym

class SolutionSystem:
    """Orchestrate a system of experiments.

    Instance Variables:
        ignore: A list of sympy symbols to ignore
    """
    def __init__(self, *args):
        """Create a new solution system given a list of list of XPSExperiments
        """
        self.systems = args

        scaling_factor, max_index = self.max_scaling_factor()
        self.scale(scaling_factor)
        self._ignore = []

        print('scaling factor:', scaling_factor, '\tmax index:', max_index)

    @property
    def ignore(self) -> List[sym.Symbol]:
        return self._ignore
    
    @ignore.setter
    def ignore(self, ignore: List[sym.Symbol]) -> None:
        self._ignore = ignore

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
