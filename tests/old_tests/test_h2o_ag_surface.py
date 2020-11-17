"""
A fast integration test for Surface CRN.

Save the resulting trajectory with month/date/year in the title into test/figures directory.
"""
import datetime
from lblcrn import *

TEST_FIGURES_DIRECTORY = "test/figures/"
date = datetime.datetime.now()
date_str = f"{date.strftime('%b').lower()}_{date.day}_{date.year}"

sm = SpeciesManager()
s = Surface("Ag",
            color="#e1e1e1",
            poscar_file="lblcrn/surface_crn/connectivity/poscar_files/POSCAR_Ni_111",
            supercell_dimensions=3,
            surface_depth=1)

h2o_star = sm.sp("H2O*", Orbital('1s', 532.2), color="orange")
h2o_multi = sm.sp("H2O_multi", Orbital('1s', 533.2), size=1, color="pink")
h2o_o_hb = sm.sp("H2O_O_hb", Orbital('1s', 531.6), color="teal")
h2o_oh_hb = sm.sp("H2O_OH_hb", Orbital('1s', 531.6), color="navy")
sm.name_be("H2O_hb", 531.6, color="blue")
o_star = sm.sp("O*", Orbital('1s', 530.0), color="red")
oh_star = sm.sp("OH*", Orbital('1s', 530.9), color="magenta")
o_star_threef = sm.sp(o_star, site=s.threefold)
oh_star_threef = sm.sp(oh_star, site=s.threefold)
o_star_twof = sm.sp(o_star, site=s.twofold)
oh_star_twof = sm.sp(oh_star, site=s.twofold)
o2_g = sm.sp("O2_gas", Orbital('1s', 535.0), color=(137, 93, 109))

rsys = RxnSystem(
    sm,
    s,

    SurfaceRevRxn([o_star], [h2o_o_hb], 1.238045, 0.127713),
    SurfaceRevRxn([oh_star], [h2o_oh_hb], 0.526204, 0.30048),
    SurfaceRevRxn([h2o_o_hb, s.threefold], [oh_star, oh_star_threef], 6.220646, 0.160755),
    #     SurfaceRevRxn([h2o_oh_hb, s.threefold], [sm.drop_marker(h2o_star, "H2O*_through_OH_hb", color="yellow"), oh_star_threef], 0.299507, 1),
    #     SurfaceRevRxn([h2o_o_hb, s.threefold], [sm.drop_marker(h2o_star, "H2O*_through_O_hb", color="#893101"), o_star_threef], 0.16713, 1),
    SurfaceRevRxn([h2o_oh_hb, s.threefold], [h2o_star, oh_star_threef], 0.299507, 1),
    SurfaceRevRxn([h2o_o_hb, s.threefold], [h2o_star, o_star_threef], 0.16713, 1),

    #     SurfaceRxn([sm.drop_marker(h2o_star, "adsorbed_H2O_gas", color="#151E3D")], [s], 0.794455),
    SurfaceRxn([h2o_star], [s], 0.794455),
    #     SurfaceRxn([s], [sm.drop_marker(h2o_star, "H2O_gas", color="#008ECC")], 0.199022),
    SurfaceRxn([s], [h2o_star], 0.199022),
    SurfaceRevRxn([h2o_oh_hb], [h2o_multi], 0.400796, 0.3945),

    SurfaceRevRxn([o_star, o_star_threef], [s, s.threefold], 0.00096, 76.964514),

    # The following two diffusion rules are useless
    # SurfaceRxn([s.threefold, o_star_threef], [o_star_threef, s.threefold], 3),
    # SurfaceRxn([s.threefold, oh_star_threef], [oh_star_threef, s.threefold], 3),

    # Revision makes use of twofold as a temporary state.
    #     SurfaceRxn([s.twofold, o_star_threef], [o_star_twof, s.threefold], 3),
    #     SurfaceRxn([s.threefold, o_star_twof], [o_star_threef, s.twofold], 1E10),
    #     SurfaceRxn([s.twofold, oh_star_threef], [oh_star_twof, s.threefold], 3),
    #     SurfaceRxn([s.threefold, oh_star_twof], [oh_star_threef, s.twofold], 1E10),

    # o_star, and oh_star diffusion accelerated to diffuse o_star_threef and oh_star_threef quickly.
    SurfaceRxn([s, o_star], [o_star, s], 20),
    SurfaceRxn([s, oh_star], [oh_star, s], 20),
    SurfaceRxn([s, h2o_star], [h2o_star, s], 3),

    SurfaceRevRxn([oh_star, s.threefold], [s, oh_star_threef], 10, 0.5),
    SurfaceRevRxn([o_star, s.threefold], [s, o_star_threef], 10, 0.5),

    Conc(o_star_threef, 36),
)

r = scrn_simulate(rsys, 20, video=False, save_trajectory=False)
r.plot_evolution(path=TEST_FIGURES_DIRECTORY,
                 title=f"trajectory_h2o_ag_{date_str}_unlimited_gas_surface",
                 save=True,
                 use_raw_data=True)
