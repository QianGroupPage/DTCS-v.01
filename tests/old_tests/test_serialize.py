"""Tests the bulk CRN

Includes test for serialization.
"""

import copy
import json
import unittest

from dtcs import *


TEST_DIR = 'TEMP' + os.sep


class TestMSON(unittest.TestCase):
    """Test MSON serialization of crn_sym objects."""

    SM_PATH = TEST_DIR + 'sm-test.json'
    RSYS_PATH = TEST_DIR + 'rsys-test.json'
    CTS_PATH = TEST_DIR + 'cts-test.json'
    XO_PATH = TEST_DIR + 'xo-test.json'
    XPS_PATH = TEST_DIR + 'xps-test.json'

    TEST_FILES = [
        SM_PATH, RSYS_PATH, CTS_PATH, XO_PATH, XPS_PATH,
    ]
    
    def setUp(self):
        # Clean up directory
        for file_path in self.TEST_FILES:
            if os.path.exists(file_path):
                os.remove(file_path)

        if not os.path.exists(TEST_DIR):
            os.mkdir(TEST_DIR)

        # Make objects to test
        sm = SpeciesManager()
        x1 = sm.sp('sinewave', [Orbital('2p 3/2', 1, 0.75),
                   Orbital('2p 1/2', 2, 0.25)])
        x2 = sm.sp('lagged', Orbital('1s', 2))
        x3 = sm.sp('laggy_chem', Orbital('1s', 2))
        x4 = sm.sp('chemical', Orbital('1s', 2))

        rsys = RxnSystem(
            sm,

            Rxn(x1, x4, k=0.01),
            Rxn(x2, x4, k=0.01),
            Rxn(x3, x4, k=0.1),

            ConcEq(x1, -sym.sin(T)),
            Conc(x2, 0.5),
            ConcDiffEq(x2, x1),
            Term(x3, x1),
        )

        cts = simulate_crn(rsys, time=5, max_step=0.1)
        xps = cts.xps

        self.sm = sm
        self.rsys = rsys
        self.cts = cts
        self.xps = xps
        
    def tearDown(self):
        # Clean up directory
        for file_path in self.TEST_FILES:
            if os.path.exists(file_path):
                os.remove(file_path)
        
        if os.path.exists(TEST_DIR):
            os.rmdir(TEST_DIR)

    def test_species_manager_serial(self):
        sm = copy.deepcopy(self.sm)

        lsm = save_and_load(sm, self.SM_PATH)

    def test_rxn_system_serial(self):
        rsys = copy.deepcopy(self.rsys)

        lrsys = save_and_load(rsys, self.RSYS_PATH)

    def test_crn_time_series_serial(self):
        cts = copy.deepcopy(self.cts)

        lcts = save_and_load(cts, self.CTS_PATH)

    def test_xps_observable_serial(self):
        xo = copy.deepcopy(self.xps.resample())

        lxo = save_and_load(xo, self.XO_PATH)

    def test_xps_experiment_serial(self):
        xps = copy.deepcopy(self.xps)

        lxps = save_and_load(xps, self.XPS_PATH)
        
    def test_bulk_crn_serial(self):
        sm = SpeciesManager()
        x1 = sm.sp('sinewave', [Orbital('2p 3/2', 1, 0.75),
                   Orbital('2p 1/2', 2, 0.25)])
        x2 = sm.sp('lagged', Orbital('1s', 2))
        x3 = sm.sp('laggy_chem', Orbital('1s', 2))
        x4 = sm.sp('chemical', Orbital('1s', 2))

        rsys = RxnSystem(
            sm,

            Rxn(x1, x4, k=0.01),
            Rxn(x2, x4, k=0.01),
            Rxn(x3, x4, k=0.1),

            ConcEq(x1, -sym.sin(T)),
            Conc(x2, 0.5),
            ConcDiffEq(x2, x1),
            Term(x3, x1),
        )

        cts = simulate_crn(rsys, time=5, max_step=0.1)
        xps = cts.xps

        with open(self.RSYS_PATH, 'w') as file:
            json.dump(rsys, file, cls=monty.json.MontyEncoder)
        with open(self.CTS_PATH, 'w') as file:
            json.dump(cts, file, cls=monty.json.MontyEncoder)
        with open(self.XPS_PATH, 'w') as file:
            json.dump(xps, file, cls=monty.json.MontyEncoder)

        with open(self.RSYS_PATH, 'r') as file:
            lrsys = json.load(file, cls=monty.json.MontyDecoder)
        with open(self.CTS_PATH, 'r') as file:
            lcts = json.load(file, cls=monty.json.MontyDecoder)
        with open(self.XPS_PATH, 'r') as file:
            lxps = json.load(file, cls=monty.json.MontyDecoder)

        # Check the right type
        self.assertEqual(type(lrsys), RxnSystem, 'Loaded type mismatch')
        self.assertEqual(type(lcts), CRNTimeSeries, 'Loaded type mismatch')
        self.assertEqual(type(lxps), XPSExperiment, 'Loaded type mismatch')

        # Sanity check
        lsm = lrsys.species_manager
        cts_two = simulate_crn(lrsys, time=5, max_step=0.1)
        self.assertEqual(sm.get('laggy_chem'), lsm.get('laggy_chem'),
                         'Loaded symbols don\'t match')
        self.assertTrue((lcts.df == cts_two.df).all().all(),
                        'Loaded a simulation mismatch.')
        self.assertTrue((lxps.df == lcts.xps.df).all().all(),
                        'Loaded a simulation mismatch.')


def save_and_load(obj, path):
    """Save the object to disk and then load it, for testing purposes."""
    with open(path, 'w') as file:
        json.dump(obj, file, cls=monty.json.MontyEncoder)
    with open(path, 'r') as file:
        return json.load(file, cls=monty.json.MontyDecoder)


if __name__ == '__main__':
    unittest.main()
