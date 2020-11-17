"""
This test case uses the exact same set of species, and rules, as originally created for the Surface CRN.
"Fe", "2F", and "3F" are taken as a species representing the proportion of empty sites on the surface.

Save the resulting trajectory with month/date/year in the title into test/figures directory.
"""

import datetime
from lblcrn.crn_sym import *
from lblcrn.surface_crn import Results
from lblcrn.experiments.solution_system import xps

TEST_FIGURES_DIRECTORY = "test/figures/"
date = datetime.datetime.now()
date_str = f"{date.strftime('%b').lower()}_{date.day}_{date.year}"

sm = SpeciesManager()
s = sm.sp("Fe", Orbital('1s', 1), color="#e1e1e1")
s_twofold = sm.sp("2F", Orbital('1s', 1), color="#e2e2e2")
s_threefold = sm.sp("3F", Orbital('1s', 1), color="#e3e3e3")
n = sm.sp("N", Orbital('1s', 530), color="lime")
n_twofold = sm.sp("N_2F", Orbital('1s', 530))
n_threefold = sm.sp("N_3F", Orbital('1s', 530))
sm.name_be("N", 530, color="lime")
n2 = sm.sp("N2", Orbital('1s', 531), color="green")
n2_threefold = sm.sp("N2_3F", Orbital('1s', 531))
sm.name_be("N2", 531, color="green")
nh = sm.sp("NH", Orbital('1s', 532), color="orange")
nh_twofold = sm.sp("NH_2F", Orbital('1s', 532))
sm.name_be("NH", 532, color="orange")
nh2 = sm.sp("NH2", Orbital('1s', 533), color="blue")
nh2_threefold = sm.sp("NH2_2F", Orbital('1s', 533))
sm.name_be("NH2", 533, color="blue")
nh3 = sm.sp("NH3", Orbital('1s', 534), color="purple")
h = sm.sp("H", Orbital('1s', 535), color="pink")
h_threefold = sm.sp("H_3F", Orbital('1s', 535))
sm.name_be("NH2", 535, color="pink")

rsys = RxnSystem(
    sm,

    # 1, 2
    RevRxn(nh2_threefold + h, s_threefold + nh3, 4.72E+02, 5.88E+04),
    # 3, 4
    RevRxn(n_twofold + h_threefold, nh_twofold + s_threefold, 1.84E+06, 7.20E+09),
    # 5, 6
    RevRxn(nh_twofold + h_threefold, s_twofold + nh2_threefold, 1.60E+11, 1.10E+06),
    # 7, 8
    RevRxn(s, n2, 7.32E+06, 4.10E+07),
    # 9, 10
    RevRxn(n2 + s_threefold, s + n2_threefold, 5.69E+10, 1.81E+09),
    # 11, 12
    RevRxn(n2_threefold + s_twofold, n_threefold + n_twofold, 2.56E+09, 2.52E+01),
    Rxn(n_threefold + s_twofold, s_threefold + n_twofold, 1E+20),
    # 13, 14
    RevRxn(s_threefold + s, h_threefold + h, 7.67E+08, 7.67E+08),
    # 15
    # Rxn(s, nh3, 6.26E+03)
    Rxn(nh3, s, 3.29E+05),
    # 16, 17
    RevRxn(nh2_threefold + s, s_threefold + nh2, 1.43E+10, 1.41E+13),
    # 18, 19
    RevRxn(nh_twofold + s, s_twofold + nh, 1.46E+07, 1.41E+13),
    # 20, 21
    RevRxn(n_twofold + s, s_twofold + n, 1.46E+07, 1.41E+13),
    # 22, 23
    RevRxn(h_threefold + s, s_threefold + h, 1.43E+10, 1.41E+13),

    Conc(n2, 5 / 48),
    Conc(nh2, 2 / 48),
    Conc(s, 1),
    Conc(s_twofold, 174/48),
    Conc(s_threefold, 127/48)
)

s, ts = xps.simulate_xps_with_cts(rsys, time=1E-8)
r = Results(None, rsys, df=ts.df.rename(columns=lambda c: str(c)))
r.plot_evolution(path=TEST_FIGURES_DIRECTORY,
                 title=f"trajectory_habor_bosch_{date_str}_unlimited_gas_bulk",
                 y_label="Coverage",
                 save=True,
                 use_raw_data=True)
