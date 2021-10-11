"""
This test case uses the exact same set of species, and rules, as originally created for the Surface CRN.
"Ag" and "3F" are taken as a species representing the proportion of empty sites on the surface.
"O2_gas" is left on the surface, so that its count or concentration can be controlled by the system.

Save the resulting trajectory with month/date/year in the title into test/figures directory.
"""
import datetime
from lblcrn.spec.crn import *
from lblcrn.sim.surface_crn import Results
from lblcrn.twin.solution_system import xps

TEST_FIGURES_DIRECTORY = "test/figures/"
date = datetime.datetime.now()
date_str = f"{date.strftime('%b').lower()}_{date.day}_{date.year}"

sm = SpeciesManager()
s = sm.sp("Ag", Orbital('1s', 1), color="#e1e1e1")
s_threefold = sm.sp("3F", Orbital('1s', 1), color="#e5e3e5")
h2o_star = sm.sp("H2O*", Orbital('1s', 532.2), color="orange")
h2o_multi = sm.sp("H2O_multi", Orbital('1s', 533.2), size=1, color="pink")
h2o_o_hb = sm.sp("H2O_O_hb", Orbital('1s', 531.6), color="teal")
h2o_oh_hb = sm.sp("H2O_OH_hb", Orbital('1s', 531.6), color="navy")
sm.name_be("H2O_hb", 531.6, color="blue")
o_star = sm.sp("O*", Orbital('1s', 530.0), color="red")
oh_star = sm.sp("OH*", Orbital('1s', 530.9), color="magenta")
o_star_threef = sm.sp("O*_3F", Orbital('1s', 530.0))
oh_star_threef = sm.sp("OH*_3F", Orbital('1s', 530.9))
o_star_twof = sm.sp("O*_2F", Orbital('1s', 530.0))
oh_star_twof = sm.sp("OH*_2F", Orbital('1s', 530.9))
sm.name_be("O*", 530.0, color="red")
sm.name_be("OH*", 530.9, color="magenta")
o2_g = sm.sp("O2_gas", Orbital('1s', 535.0), color=(137, 93, 109))

rsys = RxnSystem(
    sm,

    RevRxn(o_star, h2o_o_hb, 1.238045, 0.127713),
    RevRxn(oh_star, h2o_oh_hb, 0.526204, 0.30048),
    RevRxn(h2o_o_hb + s_threefold, oh_star + oh_star_threef, 6.220646, 0.160755),
    RevRxn(h2o_oh_hb + s_threefold, h2o_star + oh_star_threef, 0.299507, 1),
    RevRxn(h2o_o_hb + s_threefold, h2o_star + o_star_threef, 0.16713, 1),

    Rxn(h2o_star, s, 0.794455),
    Rxn(s, h2o_star, 0.199022),
    RevRxn(h2o_oh_hb, h2o_multi, 0.400796, 0.3945),

    RevRxn(o_star, o2_g, 0.00096, 76.964514),

    Rxn(s + o_star, o_star + s, 30),
    Rxn(s + oh_star, oh_star + s, 30),
    Rxn(s + h2o_star, h2o_star + s, 3),

    RevRxn(oh_star + s_threefold, s + oh_star_threef, 10, 0.5),
    RevRxn(o_star + s_threefold, s + o_star_threef, 10, 0.5),

    Conc(o_star_threef, 0.25),
    Conc(s, 1),
    Conc(s_threefold, 336/144 - 0.25),
)

s, ts = xps.simulate_xps_with_cts(rsys, time=20)
r = Results(None, rsys, df=ts.df.rename(columns=lambda c: str(c)))
r.plot_evolution(path=TEST_FIGURES_DIRECTORY,
                 title=f"trajectory_h2o_ag_{date_str}_limit_gas_bulk",
                 y_label="Coverage (1.0=144 counts)",
                 save=True,
                 use_raw_data=True)
