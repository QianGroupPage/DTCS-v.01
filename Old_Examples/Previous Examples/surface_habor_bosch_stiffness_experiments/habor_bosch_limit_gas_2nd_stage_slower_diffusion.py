"""
Habor-Bosch system for 200s, a stiff system where the simulation time would be very long.

Resulting gas from n2 and h2 desorption/adsorption are tracked in the surface system.

Assume the system has achieved relative stability after 1E-5s.  Fast diffusion rates are scaled by 4 to 6 orders of
magnitudes.
"""
from dtcs import *

sm = SpeciesManager()
s = Surface("Fe", color="#e1e1e1", poscar_file="dtcs/surface_crn/connectivity/poscar_files/POSCAR_Fe_111",
            supercell_dimensions=1,
            surface_depth=1.7)

n = sm.sp("N", Orbital('1s', 530), color="lime")
n_twofold = sm.sp(n, site=s.twofold)
n_threefold = sm.sp(n, site=s.threefold)
n2 = sm.sp("N2", Orbital('1s', 530), color="green")
n2_threefold = sm.sp(n2, site=s.threefold)
n2_g = sm.sp("N2_gas", Orbital('1s', 530), color="#043927")

nh = sm.sp("NH", Orbital('1s', 530), color="orange")
nh_twofold = sm.sp(nh, site=s.twofold)
nh2 = sm.sp("NH2", Orbital('1s', 530), color="blue")
nh2_threefold = sm.sp(nh2, site=s.threefold)
nh3 = sm.sp("NH3", Orbital('1s', 530), color="purple")

h = sm.sp("H", Orbital('1s', 530), color="pink")
h_threefold = sm.sp(h, site=s.threefold)
h2_g = sm.sp("H2_gas", Orbital('1s', 530), color="#FE7F9C")

# s_nh3_gas = measure(s, name="NH3_gas")


rsys = RxnSystem(
    sm,
    s,

    # 1, 2 Formation of adsorbed NH3
    SurfaceRevRxn([nh2_threefold, h], [s.threefold, nh3], 4.72E+02, 5.88E+04),
    # 3, 4 Formation of NH
    SurfaceRevRxn([n_twofold, h_threefold], [nh_twofold, s.threefold], 1.84E+06, 7.20E+09),
    # 5, 6
    SurfaceRevRxn([nh_twofold, h_threefold], [s.twofold, nh2_threefold], 1.60E+11, 1.10E+06),
    # 7, 8 Desorption/adsorption
    SurfaceRevRxn([n2_g], [n2], 7.32E+06, 4.10E+07),
    # 9, 10
    SurfaceRevRxn([n2, s.threefold], [s, n2_threefold], 5.69E+10, 1.81E+09),
    # 11, 12
    SurfaceRevRxn([n2_threefold, s.twofold], [n_threefold, n_twofold], 2.56E+09, 2.52E+01),
    # Changed down from 1E+20 to 1E+14 to potentially increase speed.
    SurfaceRxn([n_threefold, s.twofold], [s.threefold, n_twofold], 1E+14),
    # 13, 14 Desorption/adsorption
    SurfaceRevRxn([s.threefold, h2_g], [h_threefold, h], 7.67E+08, 7.67E+08),
    # 15
    # SurfaceRxn([s], [nh3], 6.26E+03)
    SurfaceRxn([nh3], [sm.drop_marker(s, "NH3_gas", color="red")], 3.29E+05),
    # 16, 17  diffusion rule has faster side -= 4 - 2, slower side -=2 - 2
    SurfaceRevRxn([nh2_threefold, s], [s.threefold, nh2], 1.43E+6, 1.41E+7),
    # 18, 19
    SurfaceRevRxn([nh_twofold, s], [s.twofold, nh], 1.46E+03, 1.41E+7),
    # 20, 21
    SurfaceRevRxn([n_twofold, s], [s.twofold, n], 1.46E+3, 1.41E+7),
    # 22, 23
    SurfaceRevRxn([h_threefold, s], [s.threefold, h], 1.43E+6, 1.41E+7),

    # Rates at the beginning og skipped stages.
    # Conc(n2, 5),
    # Conc(nh2, 2),
    # Conc(h2_g, 5),
    # Conc(n2_g, 5),
    Conc(n_twofold, 18),
    Conc(nh2_threefold, 6),
    Conc(h_threefold, 15),
    Conc(h, 1)
)


results = scrn_simulate(rsys, 200, video=False, spectra_in_video=True, spectra_average_duration=3, ensemble_size=1,
                        save_trajectory=True, trajectory_path="habor_bosch_limit_gas_2nd_stage_slower_diffusion", compress_trajectory=True,
                        section_length=5E5)
