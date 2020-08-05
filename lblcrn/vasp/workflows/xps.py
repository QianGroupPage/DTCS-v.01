"""
XPS Workflows

TODO(Andrew)
"""

from uuid import uuid4

from atomate.vasp.config import DB_FILE, VASP_CMD
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from atomate.vasp.powerups import add_additional_fields_to_taskdocs
from atomate.utils.utils import get_fws_and_tasks
from fireworks import Firework, Workflow
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

from lblcrn import _echo
from lblcrn.vasp.firetasks.parse_outputs import CoreEigenToDb
from lblcrn.vasp.fireworks.core import CRNSimulateFW, XPSSimulateFW

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'


def get_wf_simulate_xps(  # TODO(Andrew) Typehints
        surface_species,
        reaction_system,
        vasp_cmd=VASP_CMD,
        db_file=DB_FILE,
):
    """
    Returns an XPS simulation workflow.

    The workflow is structured:
    - (A) Fireworks for VASP calculation of core eigenenergies, two for each
            species. First is the structural relaxation/optimization, and
            second is the static binding energy calculation.
    - (B) Firework: collect all the outputted core eigenenergies into one
            central database.
    - (C) Firework: a chemical reaction network (CRN) simulation of the
            species; this generates simulated concentrations. It can be a
            bulk (homogeneous) simulation or a surface simulation.
    - (D) Firework: use the core eigenenergies and teh simulated
            concentrations to simulate an XPS experiment.
    TODO(Andrew) This method should also handle (1) multiple simulations (.e.g
     if we're using Rithvik's solution system stuff for many temperatures and
     pressures) and (2) a component Firework which waits for experimental
     data to be created and marked as 'approved', which it will then
     interpret. I'm going to KISS this for now, though.

    The dependency structure is:
    A -> B -> D <- C

    Args: TODO(Andrew)

    Returns:
        A Workflow as specified in the docstring.
    """

    # Generate workflow metadata to keep track of what is and isn't part of
    #  this workflow. This is used for database searching internally.
    wf_uuid = str(uuid4())
    wf_uuid_head = wf_uuid[:4]
    wf_name = f'Simulate XPS #{wf_uuid_head}'
    wf_meta = {
        'wf_uuid': wf_uuid,
        'wf_name': wf_name,
    }

    _echo.echo(f'Creating Workflow "{wf_name}" with uuid {wf_uuid}...')

    # Make the fireworks
    fws = []

    # For each surface specie, calculate with VASP the core eigenenergies;
    #  this adds two fireworks per surface specie, a relax and a static run.
    fws_calc_eigen = []
    for specie_index, surface_specie in enumerate(surface_species):
        structure = surface_specie.structure
        name = surface_specie.name

        # The structural relaxation firework.
        relax_vinput = MPRelaxSet(
            structure=structure,
            # user_incar_settings=...
        )
        relax_fw = OptimizeFW(
            structure=structure,
            name=f'{name} relaxation',
            vasp_input_set=relax_vinput,
            vasp_cmd=vasp_cmd,
            # override_default_vasp_params=None,
            db_file=db_file,
        )

        # The static VASP firework to calculate core eigenenergies
        core_eigen_incar_settings = {
            'ICORELEVEL': 1  # TODO: Option to make it 2
            # TODO: Update with user-supplied settings
        }
        calc_eigen_vinput = MPStaticSet(
            structure=structure,
            user_incar_settings=core_eigen_incar_settings
        )
        calc_eigen_fw = StaticFW(
            structure=structure,  # Only used for labeling
            name=f'{name} core eigenenergy calculation',
            vasp_input_set=calc_eigen_vinput,
            vasp_cmd=vasp_cmd,
            db_file=db_file,
            prev_calc_loc=True,
            parents=relax_fw,
            # TODO: If we can get atomate to add it, vasptodb_kwargs is
            #  where we'd put a '{parse_core_eigenenergies: True}' option.
            vasptodb_kwargs={
                'additional_fields': {
                    'core_eigen_calc_index': specie_index,
                },
                'extra_outcar_read': ['core_state_eigen'],
            },
        )

        # Add all the fireworks to the master list. We keep track of
        #  fws_calc_eigen separately because they are all parents of
        #  fw_aggregate_core_eigens.
        fws.append(relax_fw)
        fws.append(calc_eigen_fw)
        fws_calc_eigen.append(calc_eigen_fw)

    # Make one firework to aggregate the core eigenenergies into one place.
    fw_aggregate_core_eigens = Firework(
        tasks=CoreEigenToDb(
            calc_count=len(surface_species),
            db_file=db_file,
            wf_uuid=wf_uuid,
            wf_meta=wf_meta,
            # TODO: Perhaps add a to_db parameter?
        ),
        name=f'CoreEigenToDb for workflow {wf_name}',
        parents=fws_calc_eigen,
    )
    fws.append(fw_aggregate_core_eigens)

    # Make the firework(s) which run the simulations of the concentrations of
    #  the species. This can be a bulk (homogeneous) or surface simulation.
    fw_crn_simulate = CRNSimulateFW(
        reaction_system=reaction_system,
        sim_type='bulk',
        sim_options={'time': 10, 'max_step': 0.01},
        db_file=db_file,
        wf_meta=wf_meta,
    )
    fws.append(fw_crn_simulate)

    # Make the firework(s) which simulate XPS experiments using the simulated
    #  core eigenenergies and the simulated surface species' concentrations.
    #fw_xps_simulate = XPSSimulateFW(
    #    name=None,  # TODO
    #    parents=[fw_aggregate_core_eigens, fw_crn_simulate],
    #)
    #fws.append(fw_xps_simulate)

    # TODO(Andrew): Experiment dummy and comparison fireworks.

    # TODO: Apply common/necessary powerups, both user-supplied and built-in.
    #  Use atomate.add_common_powerups and add_additional_fields_to_taskdocs
    #  to add the metadata.

    # Make the workflow
    wf = Workflow(
        fireworks=fws,
        name=wf_name,
        metadata=wf_meta,
    )
    # Add the wf metadata to all VaspToDb's outputs so that CoreEigenToDb
    #  knows what to look for.
    wf = add_additional_fields_to_taskdocs(wf, {
        'wf_meta': wf_meta
    })

    return wf
