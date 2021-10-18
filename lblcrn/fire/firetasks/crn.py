"""
Firetasks for Chemical Reaction Network simulation

TODO(Andrew)
"""

import datetime
import json
import os

from atomate.utils.utils import env_chk, get_logger
from fireworks import FiretaskBase, FWAction, explicit_serialize
from monty.json import MontyDecoder, MontyEncoder

from lblcrn import lblcrn_echo_on
from lblcrn.twin.core import simulate

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'

_logger = get_logger(__name__)


@explicit_serialize
class BulkCRNSim(FiretaskBase):
    """
    TODO(Andrew)
    """

    required_params = ['reaction_system']
    optional_params = ['sim_options']

    def run_task(self, fw_spec):

        _logger.info(f'Running bulk CRN simulation in {os.getcwd()}')

        decode = MontyDecoder().process_decoded

        rsys = decode(self.get('reaction_system'))
        sim_options = self.get('sim_options', {})

        with lblcrn_echo_on():  # TODO: Does this even work?
            cts = simulate(
                rsys=rsys,
                **sim_options,
            )[1]

        # Make the output dict
        out_dict = {
            'cts': cts.as_dict(),
        }

        # Output MSON to cwd.
        out_path = os.path.join(os.getcwd(), 'crn_time_series.json')
        with open(out_path, 'w') as file:
            json.dump(obj=out_dict,
                      fp=file,
                      cls=MontyEncoder, )

        _logger.info(f'Saved bulk CRN simulation in {out_path}')

        return FWAction()
