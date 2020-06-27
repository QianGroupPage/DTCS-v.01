"""Tests lblcrn.crn_sym

Includes test for serialization.
"""

import json
import unittest

from lblcrn import *


TEST_DIR = 'TEMP' + os.sep


class TestMSON(unittest.TestCase):
    """Test MSON serialization of crn_sym objects."""

    RSYS_PATH = TEST_DIR + 'rsys-test.json'

    TEST_FILES = [
        RSYS_PATH,
    ]
    
    def setUp(self):
        if not os.path.exists(TEST_DIR):
            os.mkdir(TEST_DIR)
        
    def tearDown(self):
        for file_path in self.TEST_FILES:
            if os.path.exists(file_path):
                os.remove(file_path)
        
        if os.path.exists(TEST_DIR):
            os.rmdir(TEST_DIR)
        
    def test_file_io(self):
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

        with open(self.RSYS_PATH, 'w') as file:
            json.dump(rsys, file, cls=monty.json.MontyEncoder)

        with open(self.RSYS_PATH, 'r') as file:
            lrsys = json.load(file, cls=monty.json.MontyDecoder)

        # Check the right type
        self.assertEqual(type(lrsys), RxnSystem, 'Loaded type mismatch')

        # Sanity check
        lsm = lrsys.species_manager
        cts_one = simulate_crn(rsys, time_max=5, max_step=0.01)
        cts_two = simulate_crn(lrsys, time_max=5, max_step=0.01)
        self.assertEqual(sm.get('laggy_chem'), lsm.get('laggy_chem'),
                         'Loaded symbols don\'t match')
        self.assertTrue((cts_two.df == cts_one.df).all().all(),
                        'Loaded a simulation mismatch.')
    

if __name__ == '__main__':
    unittest.main()
