"""
Firetasks for parsing VASP outputs

# TODO(Andrew)
"""

import datetime
import json
import os

from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import env_chk, get_logger
from atomate.vasp.database import VaspCalcDb
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks import FiretaskBase, FWAction, explicit_serialize
from monty.json import MontyDecoder, jsanitize
import pymongo

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'

_logger = get_logger(__name__)


@explicit_serialize
class MsonToDb(FiretaskBase):  # TODO: Move this out of the vasp module.
    """
    Copy Monty JSON data from a file to a database (or to another .json file).

    Designed to be inheritable, override run_task and use ._argparse and
    ._write_out to easily make your custom database writings.

    Required Parameters:
        in_fname (str): The json file to read from.

    Optional Parameters:
        calc_dir (str): path to dir (on current filesystem) that contains
            output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name.
        db_file (str): The atomate database config file to access the database
            with. Required if you intend to write out to a database.
        db_name (str): The database collection to write out to. Required if you
            intend to write out to a database.
        to_json (bool): If True, it will write to the json file specified by
            out_dir and out_fname. Default: False
        out_dir (str): The directory to write out to if you aren't using a
            database. Defaults to current directory.
        out_fname (str): The filename to write out to if you aren't using a
            database. Defaults to in_fname.
        wf_meta (dict): Workflow metadata to include in the output.
        wf_uuid (str): If you don't give wf_meta, it will default wf_meta to
            {'wf_uuid': wf_uuid}.

    A more extensible version of atomate's JsonToDb.
    """

    required_params = ['in_fname']
    optional_params = ['calc_dir', 'calc_loc',
                       'db_file', 'db_name',
                       'to_json', 'out_dir', 'out_fname',
                       'wf_meta', 'wf_uuid']

    def _argparse(self, fw_spec):
        """
        Retrieves the following arguments, saving them in self as a dict:
            db_file: How to connect to the database (from atomate)
            db_name: Which collection to save the data under
            wf_meta: Workflow metadata to save
            in_fname: The filename to save out to
            in_fpath: The json file to load
        """

        # Get db_file
        db_file = env_chk(self.get('db_file'), fw_spec)

        # Get db_name
        db_name = self.get('db_name')

        # Get wf_meta
        wf_uuid = self.get('wf_uuid')
        wf_meta = {'wf_uuid': wf_uuid} if wf_uuid else None
        wf_meta = self.get('wf_meta', wf_meta)

        # Get in_fpath from calc_dir
        in_fname = self.get('in_fname')
        calc_dir = self.get('calc_dir', os.getcwd())
        if self.get('calc_loc'):
            calc_dir = get_calc_loc(self.get('calc_loc'),
                                    fw_spec['calc_locs'])['path']
        in_fpath = os.path.join(calc_dir, in_fname)

        # Get out_fpath and to_json
        to_json = self.get('to_json', not db_file)
        out_dir = self.get('out_dir', os.getcwd())
        out_fname = self.get('out_fname', in_fname)
        out_fpath = os.path.join(out_dir, out_fname)

        # Raise an error if they supply db_file but not db_name.
        if not to_json and not (db_name and db_file):
            raise TypeError('If outputting to database, both db_file and '
                            'db_name are required.')

        # Save to self as a dict
        self['db_file'] = db_file
        self['db_name'] = db_name
        self['wf_meta'] = wf_meta
        self['in_fpath'] = in_fpath
        self['to_json'] = to_json
        self['out_fpath'] = out_fpath

    def run_task(self, fw_spec):

        self._argparse(fw_spec)
        in_fpath = self.get('in_fpath')

        _logger.info(f'Getting data from {in_fpath}')

        # Retrieve MSON data
        with open(in_fpath, 'r') as file:
            mson_obj = json.load(file, cls=MontyDecoder)

        # Create output data
        out_dict = jsanitize(mson_obj.as_dict())
        self._write_out(out_dict)

    def _write_out(self, out_dict):

        wf_meta = self.get('wf_meta')
        db_file = self.get('db_file')
        db_name = self.get('db_name')
        to_json = self.get('to_json')
        out_fpath = self.get('out_fpath')

        out_dict.update({
            'last_updated': datetime.datetime.utcnow(),
        })
        if wf_meta:
            out_dict['wf_meta'] = wf_meta

        # Store the results in the database or a JSON file
        if to_json:
            _logger.info(f'Writing JSON to {out_fpath}.')
            with open(out_fpath, 'w') as file:
                json.dump(obj=out_dict,
                          fp=file,
                          default=DATETIME_HANDLER, )
        else:
            _logger.info(f'Uploading to database {db_name}.')
            # TODO: Should I use VaspCalcDb even though it isn't Vasp?
            mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
            mmdb.db[db_name].insert_one(out_dict)

        _logger.info('Retrieval successful.')

        return FWAction()


@explicit_serialize
class CoreEigenToDb(FiretaskBase):
    """
    Aggregate the core eigenvalues for many species (each one with its own
    VASP runs) into one database.

    In order for this task to find your tasks, you need to set
    core_eigen_calc_index in the additional_fields of VaspToDb. This lets
    CoreEigenToDb know which species is which, as often the chemical formula
    is not sufficient identification.

    Required Parameters:
        calc_count (int): The number of species you're calculating core
            eigenenergies for.
        db_file (str): Path to the db.json file which points to the database
            with your tasks collection and where you want your core_eigen
            collection.
        wf_uuid (str): The unique ID of the workflow, used to identify which
            tasks to trawl. Usually auto-generated by get_wf_...().

    Optional Parameters:
        wf_meta (str): Additional workflow metadata besides the wf_uuid.
    """

    JSON_DUMP_FILE = 'core_eigen.json'

    required_params = ['calc_count', 'db_file', 'wf_uuid']
    optional_params = ['wf_meta']

    def run_task(self, fw_spec):

        # Retrieve parameters
        calc_count = self.get('calc_count')
        db_file = env_chk(self.get('db_file'), fw_spec)
        wf_uuid = self.get('wf_uuid')
        wf_meta = self.get('wf_meta', {'wf_uuid': wf_uuid})

        # Retrieve the information about each run (in particular, core
        #  level eigenenergies) from the database.
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

        # For each core eigen calculation, create a profile of the most
        #  up-to-date run (a task) of that calculation.
        task_outputs = []
        for calc_index in range(calc_count):
            task = mmdb.collection.find_one(
                filter={
                    'wf_meta.wf_uuid': wf_uuid,
                    'core_eigen_calc_index': calc_index,
                },
                sort=[('last_updated', pymongo.DESCENDING)],
            )

            outcar = task['calcs_reversed'][0]['output']['outcar']

            task_profile = {
                # TODO: Add CONTCAR & output energy from relax run
                'chemsys': task['chemsys'],
                'core_energies': outcar['core_state_eigen'],
                'formula_pretty': task['formula_pretty'],
                'formula_reduced_abc': task['formula_reduced_abc'],
                'last_updated': task['last_updated'],
                'nelements': task['nelements'],
                'nsites': task['nsites'],
                'task_id': task['task_id'],
            }
            task_outputs.append(task_profile)

        # Make the output collection
        out_dict = {
            'last_updated': datetime.datetime.utcnow(),
            'task_outputs': task_outputs,
            'ntasks': calc_count,
            'wf_meta': wf_meta,
        }

        # Store the results in the database or a JSON file
        # TODO: Make this optional, in this setting the db_file is required.
        if not db_file:
            _logger.info(f'No database file found, attempting to dump to '
                         f'{os.path.join(os.getcwd(), self.JSON_DUMP_FILE)}.')
            with open(self.JSON_DUMP_FILE, 'w') as file:
                json.dump(obj=out_dict,
                          fp=file,
                          default=DATETIME_HANDLER, )
        else:
            mmdb.db['core_eigen'].insert_one(out_dict)

        _logger.info('Core eigenenergies successfully retrieved.')

        return FWAction()