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
import pymongo

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'

_logger = get_logger(__name__)


@explicit_serialize
class SimulateXPS(FiretaskBase):
    """
    TODO(Andrew)
    """

    required_params = []
    optional_params = []

    def run_task(self, fw_spec):

        #_logger.info(f'Running bulk XPS simulation in {os.getcwd()}')
        _logger.info(f'Pretending to run XPS simulation.')
