import unittest
import pandas as pd

from lblcrn.twin.core import simulate


class TestBulk(unittest.TestCase):
    def test_predator_prey(self):
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

        _, ts = simulate(rsys, time=45, max_step=1e-3)
        want = pd.read_csv("./tests/data/bulk_predator_prey.csv")
        for step in want.iterrows():
            got_x1 = ts.at(step[1][0])[0]
            got_x2 = ts.at(step[1][0])[1]
            want_x1 = step[1][1]
            want_x2 = step[1][2]

            x1_correct, x2_correct = False, False
            if got_x1 <= want_x1*1.003 and got_x1 >= want_x1*0.997:
                x1_correct = True
            if got_x2 <= want_x2*1.003 and got_x2 >= want_x2*0.997:
                x2_correct = True
            self.assertTrue(x1_correct, f"Expected {want_x1}, got {got_x1}")
            self.assertTrue(x2_correct, f"Expected {want_x2}, got {got_x2}")

    def test_ag(self):
        sm = SpeciesManager()

        y1 = sm.sp('H2Og', Orbital('1s', 535.0))
        x2 = sm.sp('H2O*', Orbital('1s', 532.2))
        x3 = sm.sp('OH*', Orbital('1s', 530.9))
        x4 = sm.sp('O*', Orbital('1s', 530.0))
        x53 = sm.sp('OH.H2O_hb', Orbital('1s', 531.6))
        x54 = sm.sp('O.H2O_hb', Orbital('1s', 531.6))
        x6 = sm.sp('multiH2O', Orbital('1s', 533.2))
        x7 = sm.sp('O2g', Orbital('1s', 535.0))

        rsys = RxnSystem(
            Rxn(x4 + y1, x54, 3.207654),
            Rxn(x3 + y1, x53, 1.363342),
            RevRxn(x54, x3 + x3, 6.220646,0.160755),
            Rxn(x53, x2 + x3, 0.299507),
            Rxn(x54, x2 + x4, 0.167130),
            Rxn(x2, y1, 1.939313),
            Rxn(y1, x2, 0.515646),
            Rxn(x53, y1 + x3, 0.733491),
            Rxn(x54, x4 + y1, 0.311754),
            Rxn(x53 + y1, x6, 1.038423),
            Rxn(x6, x53 + y1, 0.962999),
            RevRxn(x4 + x4, x7, 0.002342,426.922895),
            Conc(y1,1),
            Conc(x4,0.25),
            sm
        )

        _, ts = simulate(rsys, 500, max_step=1)
        # want = pd.read_csv("./tests/data/bulk_ag.csv")
        # for step in want.iterrows():
            # got = ts.at(step[0])

            # for got_species in got:
                # correct = False
                # if got_species <= want_species*1.003 and got_species >= want_species*0.997:
                    # correct = True

                # self.assertTrue(x1_correct, f"Expected {want_x1}, got {got_x1}")

    def test_nanomolar(self):
        sm = SpeciesManager()
        a, b, c = sm.sp('a'), sm.sp('b'), sm.sp('c')

        rsys = RxnSystem(
            RevRxn(a+b, c, 1.5e6, 10),
            Rxn(2*a,b, 1e5),
            Conc(a, 1e-9),
            Conc(b, 1e-9),
            Conc(c, 1e-9),
            sm
        )
        _, ts = simulate(rsys, 3600 * 10, max_step=1)
        want = pd.read_csv("./tests/data/bulk_nanomolar.csv")

        for step in want.iterrows():
            got_x1 = ts.at(step[1][0])[0]
            got_x2 = ts.at(step[1][0])[1]
            want_x1 = step[1][1]
            want_x2 = step[1][2]

            x1_correct, x2_correct = False, False
            if got_x1 <= want_x1*1.04 and got_x1 >= want_x1*0.9:
                x1_correct = True
            if got_x2 <= want_x2*1.04 and got_x2 >= want_x2*0.9:
                x2_correct = True
            self.assertTrue(x1_correct, f"Expected {want_x1}, got {got_x1}")
            self.assertTrue(x2_correct, f"Expected {want_x2}, got {got_x2}")
