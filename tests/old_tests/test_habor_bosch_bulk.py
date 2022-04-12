"""
Test Habor Bosch system using all stiff-oriented methods.

Credits: Habor Bosch originally designed by Rithvik. Testing script by Ye.
"""

from dtcs.twin import XPSInitializationData
from dtcs.sim.surface_crn import Results
import datetime

TEST_FIGURES_DIRECTORY = "test/figures/"
date = datetime.datetime.now()
date_str = f"{date.strftime('%b').lower()}_{date.day}_{date.year}"

sm = SpeciesManager()

v_3n_nh2 = sm.sp('v_3n_nh2', Orbital('1s', 535.1))
v_h2 = sm.sp('v_h2', Orbital('1s', 532.2))
v_3n_nh2_2h = sm.sp('v_3n_nh2_2h', Orbital('1s', 530.9))
v_3n_nh3_h = sm.sp('v_3n_nh3_h', Orbital('1s', 530.2))
v_3n_h = sm.sp('v_3n_h', Orbital('1s', 530.3))
v_nh3 = sm.sp('v_nh3', Orbital('1s', 530.4), color="black")
v_2n_nh = sm.sp('v_2n_nh', Orbital('1s', 530.5))
v_2n_nh_2h = sm.sp('v_2n_nh_2h', Orbital('1s', 530.6))
v_2n_nh2_h = sm.sp('v_2n_nh2_h', Orbital('1s', 530.7))
v_2n_nh3 = sm.sp('v_2n_nh3', Orbital('1s', 534.8))
v_2n = sm.sp('v_2n', Orbital('1s', 530.9))
v_n2 = sm.sp('v_n2', Orbital('1s', 531), color="green")
v_2n_n2 = sm.sp('v_2n_n2', Orbital('1s', 531.1))
v_4n = sm.sp('v_4n', Orbital('1s', 531.2))
v_3n_nh_h = sm.sp('v_3nh_h', Orbital('1s', 531.3))
v_3n_nh2 = sm.sp('v_3nh2', Orbital('1s', 531.5))
v_2n_2h = sm.sp('v_2n_2h', Orbital('1s', 531.6))


constants = [7.67e8, 1.41e13, 2.6e6, 2.48e4, 1.03e7, 6.66e2, 1.84e6, 7.2e9, 7.67e8, 3.04e9, 1.6e11, 1.1e6, 1.14e0, 7.79e5,
            3.91e5, 3.14e3, 8.69e6, 4.1e7, 7.2e9, 1.35, 2.77e5, 6.45e8, 1.41e13, 7.67e8, 7.67e8, 1.48e4]
# constants = [10e-9 * c for c in constants]
multipliers = [1]

main = XPSInitializationData(
                'High P, High T',
                0,
                0,
                constants=constants
            )

init_data = [main]

def rsys_generator(scaled):
    rsys = RxnSystem(
        RevRxn(v_3n_nh2 + v_h2, v_3n_nh2_2h, scaled[0], scaled[1]),
        RevRxn(v_3n_nh2_2h, v_3n_nh3_h, scaled[2], scaled[3]),
        RevRxn(v_3n_nh3_h, v_3n_h + v_nh3, scaled[4], scaled[5]),
        RevRxn(v_3n_h, v_2n_nh, scaled[6], scaled[7]),
        RevRxn(v_2n_nh + v_h2, v_2n_nh_2h, scaled[8], scaled[9]),
        RevRxn(v_2n_nh_2h, v_2n_nh2_h, scaled[10], scaled[11]),
        RevRxn(v_2n_nh2_h, v_2n_nh3, scaled[12], scaled[13]),
        RevRxn(v_2n_nh3, v_2n + v_nh3, scaled[14], scaled[15]),
        RevRxn(v_2n + v_n2, v_2n_n2, scaled[16], scaled[17]),
        RevRxn(v_2n_n2, v_4n, scaled[18], scaled[19]),
        RevRxn(v_4n + v_h2, v_3n_nh_h, scaled[20], scaled[21]),
        RevRxn(v_3n_nh_h, v_3n_nh2, scaled[22], scaled[23]),
        RevRxn(v_2n + v_h2, v_2n_2h, scaled[24], scaled[25]),

        Conc(v_n2, 5),
        Conc(v_h2, 15),
        Conc(v_3n_nh2, 1), # 2-20 range
        sm
    )
    return rsys


methods = ['Radau', 'BDF', 'LSODA']

for method in methods:
    title = f"habor_bosch_bulk_crn_trajectory_{method.lower()}_method"
    rsys = rsys_generator(constants)
    s, ts = xps.simulate_xps_with_cts(rsys, method=method, time=200, title="Habor Bosch")

    r = Results(None, rsys, df=ts.df.rename(columns=lambda c: str(c)))

    start = ts.at(0)
    names_in_figure = []
    for i in range(len(start)):
        if str(start.keys()[i]) not in ["v_h2", "v_n2", "v_3nh2"]:
            names_in_figure.append(str(start.keys()[i]))

    r.plot_evolution(names_in_figure=names_in_figure,
                     end_time=200,
                     title=title,
                     save=True,
                     use_raw_data=True,
                     y_label="Coverage")


