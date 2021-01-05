"""
Docstring
"""
from lblcrn.crn_sym import *
from lblcrn.experiments.simulate import simulate
from lblcrn.experiments.xps_io import read_exp

sm = SpeciesManager()
x1 = sm.sp('x')
x2 = sm.sp('y')

rsys = RxnSystem(
    Rxn(x1 + x2, 2 * x2, 1.5),
    Rxn(x1, 2 * x1, 1),
    Rxn(x2, None, 1),
    Conc(x1, 2),
    Conc(x2, 1),
    sm
)

_, ts = simulate(rsys, time=40, max_step=1e-2)
ts.plot()
