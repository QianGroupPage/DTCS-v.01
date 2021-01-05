"""
HIV Bulk CRN Time Series
"""
from lblcrn.crn_sym import *
from lblcrn.experiments.simulate import simulate
from lblcrn.experiments.xps_io import read_exp

sm = SpeciesManager()
v = sm.sp('virus')
h = sm.sp('healthy')
inf = sm.sp('infected')

rsys = RxnSystem(
    sm,
    Rxn(h + v, inf , 2e-7), 
    Rxn(inf, None, 0.5),
    Rxn(h, None, 0.2),
    Rxn(v, None, 5),
    Rxn(None, h, 1e5),
    Rxn(inf, v + inf, 100),
    
    Conc(h, 1000000),
    Conc(v, 100),
)

xps, ts = simulate(rsys, 45, max_step=0.01)

ts.plot()

import matplotlib.pyplot as plt

plt.yscale("log")
plt.xlim(0,40)
