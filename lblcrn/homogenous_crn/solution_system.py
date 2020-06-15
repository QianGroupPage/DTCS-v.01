import numpy as np
import matplotlib.pyplot as plt

class SolutionSystem:
    def __init__(self, *args):
        self.systems = args

    def process(self, ignore=[]):
        for sys in self.systems:
            for s in sys: 
                s.ignore = ignore
                s.process(scale=False)

        scaling_factor, max_index = self.max_scaling_factor()
        print('scaling factor:', scaling_factor, '\tmax index:', max_index)
        self.scale(scaling_factor)

    def max_scaling_factor(self) -> int:
        """Find the maximum peak and return the scaling factor for that data and the system index.
        """
        max_env = max(self.systems[0][0].envelope)
        max_index = 0
        max_exp = max(self.systems[0][0].xps.intensity)

        # Find max peak of systems
        for i in range(0, len(self.systems)):
            for j in range(0, len(self.systems[i])):
                potential = max(self.systems[i][j].xps.intensity)
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
                new_envelope = []
                new_dists = []

                for v in s.envelope:
                    new_envelope.append(v * scaling_factor)

                for d in s.distributions:
                    new_dists.append(d * scaling_factor)

                s.envelope = new_envelope
                s.distributions = new_dists
                s.resampled_envelope = np.array(s.envelope)
    
    def __getitem__(self, index):
        return self.systems[index]
    
    def plot(self, index, rows=0, cols=0):
        fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(60,60))
        for i in range(len(self.systems[index])):
            self.systems[index][i].plot_gaussian(envelope=True, resample_envelope=True, overlay=True, ax=axes[int(i/5), int(i%5)], title=('Eq: ' + str(int(i/5)) + ' Const: ' + str(i % 5)))
        plt.show()
