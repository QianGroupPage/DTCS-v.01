import argparse
import sys
import os

from lblcrn.crn_sym import *
from lblcrn.experiments.simulate import simulate
from lblcrn.experiments.xps_io import read_exp
from lblcrn.experiments.storage import CRNStorage

sm = SpeciesManager()

y1 = sm.sp('H2Og', Orbital('1s', 535.0))
x2 = sm.sp('H2O', Orbital('1s', 532.2))
x3 = sm.sp('OH', Orbital('1s', 530.9))
x4 = sm.sp('O', Orbital('1s', 530.0))
x53 = sm.sp('OH.H2O_hb', Orbital('1s', 531.6))
x54 = sm.sp('O.H2O_hb', Orbital('1s', 531.6))
x6 = sm.sp('multiH2O', Orbital('1s', 533.2))
x7 = sm.sp('O2g', Orbital('1s', 535.0))


constants = [3.207654,1.363342,6.220646,0.160755,0.299507,0.167130,1.939313,0.515646,0.733491,0.311754,1.038423, 0.962999,0.002342,426.992895]
multipliers = [1]

def rsys_generator(scaled):
    rsys = RxnSystem(
        Rxn(x4 + y1, x54, scaled[0]),
        Rxn(x3 + y1, x53, scaled[1]),
        Rxn(x54, x3 + x3, scaled[2]),
        Rxn(x3 + x3, x54, scaled[3]),
        Rxn(x53, x2 + x3, scaled[4]),
        Rxn(x54, x2 + x4, scaled[5]),
        Rxn(x2, y1, scaled[6]),
        Rxn(y1, x2, scaled[7]),
        Rxn(x53, y1 + x3, scaled[8]),
        Rxn(x54, x4 + y1, scaled[9]),
        Rxn(x53 + y1, x6, scaled[10]),
        Rxn(x6, x53 + y1, scaled[11]),
        Rxn(x4 + x4, x7, scaled[12]),
        Rxn(x7, x4 + x4, scaled[13]),
        Conc(y1, 1),
        Conc(x4, 0.25),
        sm
    )
    return rsys

if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument("--equation", "-e", type=int, default=0, help="Specify the index of which equation to modify")
    # args = parser.parse_args()
    equation = int(sys.argv[1])
    cs = CRNStorage(uri=os.getenv("MONGO_URI"))

    for i in range(len(multipliers)):
        scaled = list(constants)
        scaled[equation] *= multipliers[i]

        rsys = rsys_generator(scaled)
        
        xps, ts = simulate(rsys, time=500, title="Eq: " + str(equation) + "Constant: " + str(i))

        # Check if the system has already been simulated
        existing = cs.load_from_rsys_id(rsys)
        if len(existing) > 0:
            continue

        try:
            cs.store(xps, ts)
            print('Solved for ('+str(equation)+', '+str(i)+')')
        except Exception as e:
            print("unable to solve for" + str(equation) + ", " + str(i))
