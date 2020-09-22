"""
XPS Workflows

TODO(Andrew)
"""

from uuid import uuid4

from atomate.vasp.config import DB_FILE, VASP_CMD
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from atomate.vasp.powerups import add_additional_fields_to_taskdocs
from fireworks import Firework, Workflow
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

from lblcrn import _echo
from lblcrn.fire.firetasks.glue_tasks import ForwardCoreEigen, ForwardSimConcs
from lblcrn.fire.fireworks.xps import XPSSimulateFW
from lblcrn.fire.fireworks.crn import CRNSimulateFW

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'


def get_wf_simulate_xps(  # TODO(Andrew) Typehints
        surface_species,
        xps_spec=None,
        crn_spec=None,
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

    Specification Schema:
        surface_species: a list of dicts as with keys:
            name (str): Used for identification.
            orbitals (dict): A dict of '<orbital name>': <binding energy>. Either
                give this or give a structure and incar to have VASP generate it.
            structure: TODO
            incar: TODO
        xps_spec:
            TODO: see XPSExperiment
        crn_spec:
            TODO

    The dependency structure is:
    A -> B -> D <- C

    Args: TODO(Andrew)

    Returns:
        A Workflow as specified in the docstring.
    """

    # Sanitize inputs
    surface_species = {specie['name']: specie for specie in surface_species}
    xps_spec = xps_spec or {}
    xps_spec['species'] = xps_spec.get('species', {})
    xps_spec['species'].update(surface_species)
    crn_spec = crn_spec or {}
    #if not type(crn_spec) == dict:
    #    crn_spec = jsanitize(crn_spec.as_dict())

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
    vasp_species = []  # TODO: Use this for ID
    calc_eigen_task_labels = {}
    for specie_index, (name, surface_specie) in enumerate(xps_spec['species'].items()):  # TODO: Stop using specie_index
        structure = surface_specie.get('structure', None)

        # If this species doesn't have a structure to do VASP on but provides
        #  orbitals manually, continue. However, if they provide neither,
        #  raise an error.
        if structure is None:
            if surface_specie.get('orbitals', None):
                continue
            else:
                raise ValueError(f'For species {name}, neither structure nor '
                                 f'orbitals were specified.')

        # The structural relaxation firework.
        relax_vinput = MPRelaxSet(
            structure=structure,
            # user_incar_settings=...
        )
        relax_fw = OptimizeFW(
            structure=structure,
            name=f'relaxation',
            vasp_input_set=relax_vinput,
            vasp_cmd=vasp_cmd,
            # override_default_vasp_params=None,
            db_file=db_file,
        )

        # The static VASP firework to calculate core eigenenergies
        calc_eigen_task_label = f'{name} core eigenenergy calculation'
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
            name=calc_eigen_task_label,  # Used for identification
            # TODO: vasp_input_set won't get used if we define parents,
            #  it seems that we have to use vasp_input_set_params instead.
            vasp_input_set=calc_eigen_vinput,
            vasp_cmd=vasp_cmd,
            db_file=db_file,
            prev_calc_loc=True,
            parents=relax_fw,
            # TODO: If we can get atomate to add it, vasptodb_kwargs is
            #  where we'd put a '{parse_core_eigenenergies: True}' option.
            vasptodb_kwargs={
                'additional_fields': {
                    'core_eigen_calc_index': specie_index,  # TODO: Stop using
                    # TODO: Instead, use wf_meta with the species name
                    #  I've decided to rely on species name as unique within
                    #  a workflow, as it needs to be unique for RxnSystem
                    #  to work anyway.
                },
                'extra_outcar_read': ['core_state_eigen'],  # TODO: This is a bodge
            },
        )

        # Add all the fireworks to the master list. We keep track of
        #  fws_calc_eigen separately because they are all parents of
        #  fw_aggregate_core_eigens.
        fws.append(relax_fw)
        fws.append(calc_eigen_fw)

        # Keep track of which surface species are using VASP
        vasp_species.append(name)
        fws_calc_eigen.append(calc_eigen_fw)
        calc_eigen_task_labels[name] = calc_eigen_task_label

    # Make one firework to aggregate the core eigenenergies into one place.
    fw_aggregate_core_eigens = Firework(
        tasks=ForwardCoreEigen(
            task_labels=calc_eigen_task_labels,
            db_file=db_file,
            wf_uuid=wf_uuid,
        ),
        name=f'ForwardCoreEigen for workflow {wf_name}',
        parents=fws_calc_eigen,
    )
    fws.append(fw_aggregate_core_eigens)

    # Make the firework(s) which run the simulations of the concentrations of
    #  the species. This can be a bulk (homogeneous) or surface simulation.
    fw_crn_simulate = CRNSimulateFW(
        reaction_system=crn_spec['reaction_network'],
        sim_type='bulk',
        sim_options={'time': 10, 'max_step': 0.01},
        db_file=db_file,
        wf_meta=wf_meta,
    )
    fws.append(fw_crn_simulate)

    # Make a firework to send the simulated concentrations from the crn
    #  simulation to XPS's spec.
    fw_aggregate_sim_concs = Firework(
        tasks=ForwardSimConcs(
            db_file=db_file,
            wf_uuid=wf_uuid,
        ),
        name=f'ForwardSimConcs for workflow {wf_name}',
        parents=fw_crn_simulate,
    )
    fws.append(fw_aggregate_sim_concs)

    # Make the firework(s) which simulate XPS experiments using the simulated
    #  core eigenenergies and the simulated surface species' concentrations.
    fw_xps_simulate = XPSSimulateFW(
        # TODO: How should this know what to do?
        spec=xps_spec,
        db_file=db_file,
        wf_meta=wf_meta,
        parents=[fw_aggregate_core_eigens, fw_aggregate_sim_concs],
    )
    fws.append(fw_xps_simulate)

    # TODO(Andrew): Experiment dummy and comparison fireworks.

    # Make the workflow
    wf = Workflow(
        fireworks=fws,
        name=wf_name,
        metadata=wf_meta,
    )

    # TODO: Apply common/necessary powerups, both user-supplied and built-in.
    #  Use atomate.add_common_powerups.

    # Add the wf metadata to all VaspToDb's outputs so that ForwardCoreEigen
    #  knows what to look for.
    wf = add_additional_fields_to_taskdocs(wf, {
        'wf_meta': wf_meta
    })

    return wf
