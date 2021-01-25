import argparse

from lblcrn.crn_sym import *
from lblcrn.experiments.simulate import simulate
from lblcrn.experiments.xps_io import read_exp

sm = SpeciesManager()

# FIXME: Declare species

constants = [] # FIXME: Set constants
multipliers = [] # FIXME: Set multipliers

def rsys_generator(scaled):
    rsys = RxnSystem(
        # FIXME: Add reactions and concentration info
        sm
    )
    return rsysh

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--equation", "-e", type=int, default=1, help="Specify the index of which equation to modify")
    args = parser.parse_args()

    sols = []
    cts = []
    sol_mult_1 = None
    cts_mult_1 = None
    for i in range(len(multipliers)):
        if multipliers[i] == 1 and sol_mult_1:
            sols.append(sol_mult_1)
            cts.append(cts_mult_1)
        else:
            scaled = list(constants)
            scaled[args.equation] *= self.multipliers[i]

            rsys = self.rsys_generator(scaled)
            xps, ts = simulate(rsys, time=self.time, title=data.title + " Eq: " + str(args.equation) + "Constant: " + str(i))

            cts.append(ts)
            sols.append(xps)
            if self.multipliers[i] == 1:
                sol_mult_1 = xps
                cts_mult_1 = ts

        print('Solved for ('+str(args.equation)+', '+str(i)+')')
        print(scaled)
        print('\n')

    # TODO: How do we want to output this data (probably jsonify it and write to a file).
