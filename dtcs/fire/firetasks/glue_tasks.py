"""
Firetasks for connecting different simulations and otherwise getting fireowrks
to work with one another.

# TODO(Andrew)
"""

from atomate.utils.utils import get_logger
from atomate.utils.utils import env_chk
from atomate.vasp.database import VaspCalcDb
from fireworks import explicit_serialize, FiretaskBase, FWAction
from monty.json import MontyDecoder
import pymongo

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'

_logger = get_logger(__name__)


@explicit_serialize
class ForwardCoreEigen(FiretaskBase):  # TODO: Rename to PushCoreEigen?
    """
    Aggregate the core eigenvalues for many species (each one with its own
    VASP runs) and forwards them to the child fireworks' specs.

    In order for this task to find your tasks, you need to set
    core_eigen_calc_index in the additional_fields of VaspToDb. This lets
    ForwardCoreEigen know which species is which, as often the chemical formula
    is not sufficient identification.

    Required Parameters:
        db_file (str): Path to the db.json file which points to the database
            with your tasks collection and where you want your core_eigen
            collection.
        task_labels (dict): A dictionary {species name: task label} used to
            identify which tasks correspond to which species' core eigenenergy
            calculations.
        wf_uuid (str): The unique ID of the workflow, used to identify which
            tasks to trawl. Usually auto-generated by get_wf_...().
    """

    required_params = ['db_file', 'task_labels', 'wf_uuid']
    # TODO: Make it able to read from JSON too, from calc_locs
    #  The previous fireworks are already passing calc_locs inot the spec,
    #  I'm just not using them for anything.
    optional_params = []

    def run_task(self, fw_spec):

        # Retrieve parameters
        db_file = env_chk(self.get('db_file'), fw_spec)
        task_labels = self.get('task_labels')
        wf_uuid = self.get('wf_uuid')

        # Retrieve the information about each run (in particular, core
        #  level eigenenergies) from the database.
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

        # For each core eigen calculation, create a profile of the most
        #  up-to-date run (a task) of that calculation.
        species_spec_update = {}
        for species_name, task_label in task_labels.items():
            task = mmdb.collection.find_one(
                filter={
                    'wf_meta.wf_uuid': wf_uuid,
                    'task_label': task_label,
                },
                sort=[('last_updated', pymongo.DESCENDING)],
            )

            species_spec = self._task_to_species_spec(task)
            species_spec_update[species_name] = species_spec

        # Make the spec update
        spec_update = {
            'species': species_spec_update,
        }

        mod_spec = _make_mod_spec(spec_update)
        return FWAction(mod_spec=[mod_spec])

    @staticmethod
    def _task_to_species_spec(task):
        """Take a VASP task's output and convert it to species spec format."""

        outcar = task['calcs_reversed'][0]['output']['outcar']
        core_eigen = outcar['core_state_eigen']
        relaxed_structure = task['output']['structure']
        sites = relaxed_structure['sites']

        # Collect the orbital information
        orbitals = []
        for (site_num, site), core_levels in zip(enumerate(sites), core_eigen):
            for orbital, energy_list in core_levels.items():
                # TODO: Convert 2p -> 2p 1/2 and
                # TODO: This will run an error with non-element sites:
                element = site['species'][0]['element']
                is_surface = element == 'Ag'  # TODO: DEBUG
                if is_surface: continue  # TODO: DEBUG

                orbitals.append({
                    'name': f'{element} {orbital} #{site_num}',
                    'element': element,
                    'orbital': orbital,
                    'binding_energy': -1 * energy_list[0],
                    'splitting': 1,  # TODO
                    # TODO: cross section
                    'site_num': site_num,
                    'is_surface': is_surface,
                })

        return {
            'orbitals': orbitals,
            # TODO: Should this go to the database or to the spec?:
            # 'relaxed_structure': relaxed_structure,
        }


@explicit_serialize
class ForwardSimConcs(FiretaskBase):
    """
    TODO
    """

    required_params = ['db_file', 'wf_uuid']
    optional_params = []

    def run_task(self, fw_spec):

        # Retrieve parameters
        db_file = env_chk(self.get('db_file'), fw_spec)
        wf_uuid = self.get('wf_uuid')

        # Access the database
        # TODO: Shouldn't I do something like CalcDb(collection='crn_sims')?
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        collection = mmdb.db['crn_sims']

        # Retrieve this workflow's most recent CRN simulation.
        cts_mson = collection.find_one(
            filter={'wf_meta.wf_uuid': wf_uuid, },
            sort=[('last_updated', pymongo.DESCENDING)],
        )

        # TODO: This does not work, look at the format in the DB for more
        #  info.
        cts = MontyDecoder().process_decoded(cts_mson['cts'])
        sim_concs = {}
        for symbol, conc in cts.at(t=-1).items():
            sim_concs[symbol.name] = conc

        # Make the spec update
        spec_update = {'sim_concs': sim_concs}

        mod_spec = _make_mod_spec(spec_update)
        return FWAction(mod_spec=[mod_spec])


def _make_mod_spec(updates: dict) -> dict:
    """
    Convert a dict into Firework's DictMod language in a sensible manner.
    """
    to_set = {}
    to_push_all = {}

    for key, value in updates.items():
        if isinstance(value, dict):
            recursed = _make_mod_spec(value)
            for key2, value2 in recursed['_set'].items():
                to_set[key + '->' + key2] = value2
            for key2, value2 in recursed['_push_all'].items():
                to_push_all[key + '->' + key2] = value2
        elif isinstance(value, list):
            to_push_all[key] = value
        else:
            to_set[key] = value
    return {
        '_set': to_set,
        '_push_all': to_push_all,
    }