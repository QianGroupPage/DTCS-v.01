"""
Firetasks for Chemical Reaction Network simulation

TODO(Andrew)
"""

import datetime
import json
import os

from atomate.utils.utils import env_chk, get_logger
from atomate.vasp.database import VaspCalcDb
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks import FiretaskBase, FWAction, explicit_serialize
from monty.json import MontyEncoder
import numpy as np
import pymongo
from sympy import Symbol

from lblcrn import lblcrn_echo_on
from lblcrn import SpeciesManager, Orbital, XPSExperiment

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'

_logger = get_logger(__name__)


@explicit_serialize
class SimulateXPS(FiretaskBase):
    """
    TODO(Andrew)
    """

    # TODO: you should be able to specify parameters here too instead of just
    #  in the spec. Or is that already possible?
    # TODO: It kind of is, task parameters are in _tasks.
    required_params = []
    optional_params = ['species', 'sim_concs']

    def run_task(self, fw_spec):

        _logger.info(f'Running bulk XPS simulation in {os.getcwd()}')

        # Retrieve simulation spec from the fw spec
        # It prioritizes inputs explicitly passed to the task (under _tasks)
        #  over inputs in the greater fw_spec.
        species_raw = self.get('species', fw_spec['species'])
        sim_concs_raw = self.get('sim_concs', fw_spec['sim_concs'])

        # TODO: This part is a bodge as there's a formatting discrepancy
        # species
        sm = SpeciesManager()
        for name, species in species_raw.items():
            orbitals = []
            for orbital in species['orbitals']:
                orbitals.append(Orbital(
                    name=orbital['name'],
                    binding_energy=orbital['binding_energy'],
                    splitting=orbital.get('splitting', 1.0)
                ))

            sm.sp(
                name=name,
                orbitals=orbitals,
                color=species.get('color', None)
            )

        # sim_concs
        sim_concs = {}
        for name, conc in sim_concs_raw.items():
            sim_concs[Symbol(name)] = conc
        # TODO: End bodge.

        # Run the XPS simulation
        with lblcrn_echo_on():  # TODO: Does this even work?
            xps = XPSExperiment(
                species_manager=sm,
                sim_concs=sim_concs,
                autoresample=False,
            )
            xps_observable = xps.resample(
                x_range=np.arange(502.5, 512.5, 0.01),  # TODO: DEBUG
                overwrite=False,
            )

        out_dict = {
            'xps_spec': xps.as_dict(),
            'xps_observable': xps_observable.as_dict(),
        }

        # Output MSON to cwd.
        out_path = os.path.join(os.getcwd(), 'xps.json')
        with open(out_path, 'w') as file:
            json.dump(obj=out_dict,
                      fp=file,
                      cls=MontyEncoder, )

        _logger.info(f'Saved XPS simulation in {out_path}')

        return FWAction()
