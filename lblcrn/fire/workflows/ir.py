"""
IR Workflows

TODO(Andrew)
"""

from uuid import uuid4

from atomate.common.firetasks.glue_tasks import CopyFilesFromCalcLoc
from atomate.vasp.config import DB_FILE, VASP_CMD
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from atomate.vasp.powerups import add_additional_fields_to_taskdocs
from fireworks import Firework, Workflow
from fireworks.user_objects.firetasks.script_task import ScriptTask, PyTask
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

from lblcrn import _echo
from lblcrn.fire.firetasks.glue_tasks import ForwardCoreEigen, ForwardSimConcs
from lblcrn.fire.fireworks.xps import XPSSimulateFW
from lblcrn.fire.fireworks.crn import CRNSimulateFW

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'


def get_wf_simulate_ir(
    species,
    vasp_cmd=VASP_CMD,
    db_file=DB_FILE,
):
    """TODO

    Args:
        species:
        vasp_cmd:
        db_file:

    Returns:
        TODO
    """
    # Sanitize input


    # Generate workflow metadata to keep track of what is and isn't part of
    #  this workflow. This is used for database searching internally.
    wf_uuid = str(uuid4())
    wf_uuid_head = wf_uuid[:4]
    wf_name = f'Simulate IR #{wf_uuid_head}'
    wf_meta = {
        'wf_uuid': wf_uuid,
        'wf_name': wf_name,
    }

    _echo.echo(f'Creating Workflow "{wf_name}" with uuid {wf_uuid}...')


def get_fws_ir_species(
    structure,
    relax_vis,
    born_vis,
    fconsts_vis,
    supercell_matrix=(1, 1, 1),
    vasp_cmd=VASP_CMD,
    db_file=DB_FILE,
    tracking_id='',
):
    """TODO

    Args:
        structure:
        relax_vis:
        born_vis:
        fconsts_vis:
        supercell_matrix:
        vasp_cmd:
        db_file:
        tracking_id:

    Returns:
        TODO
    """
    # Make the tracking ID easy to read
    if tracking_id:
        tracking_id = ' #' + tracking_id

    # Firework: VASP structural relaxation/optimization
    fv_vasp_relax_name = f'Structural Optimization' + tracking_id
    fw_vasp_relax = OptimizeFW(
        structure=structure,
        name=fv_vasp_relax_name,
        vasp_input_set=relax_vis,
        vasp_cmd=vasp_cmd,
        db_file=db_file,
    )

    # Firework: born charge VASP calculation
    fw_vasp_born_name = f'Born-Effective Charges Calc.' + tracking_id
    fw_vasp_born = StaticFW(
        structure=structure,
        name=fw_vasp_born_name,
        vasp_input_set=born_vis,
        prev_calc_loc=fv_vasp_relax_name,
        vasp_cmd=vasp_cmd,
        db_file=db_file,
        parents=fw_vasp_relax,
    )

    # Firework: force constants VASP calculation
    fw_vasp_force_consts_name = f'Force Constants Calc.' + tracking_id
    fw_vasp_force_consts = StaticFW(
        structure=structure * supercell_matrix,  # Optionally use a supercell
        name=fw_vasp_force_consts_name,
        prev_calc_loc=False,
        vasp_input_set=fconsts_vis,
        vasp_cmd=vasp_cmd,
        db_file=db_file,
    )

    # Firework: create BORN file
    fw_outcar_to_born_name = 'Create BORN file' + tracking_id
    fw_outcar_to_born = Firework(
        tasks=[
            CopyFilesFromCalcLoc(
                calc_loc=fw_vasp_born_name,
                filenames=['OUTCAR']
            ),
            ScriptTask(script='outcar-born')],
        name=fw_outcar_to_born_name,
        parents=fw_vasp_born,
    )

    # Firework: Make mask.yaml
    fw_make_mask_name = 'Make mask.yaml' + tracking_id
    fw_make_mask = Firework(
        tasks=[
            CopyFilesFromCalcLoc(
                calc_loc=fw_vasp_force_consts_name,
                filenames=['vasprun.xml']
            ),
            ScriptTask(script='phonopy -fc --mask vasprun.xml')],
        name=fw_make_mask_name,
        parents=[fw_vasp_force_consts],
    )

    # Firework: Calculate IR Spectrum
    fw_ir_spectrum = Firework(
        tasks=[
            CopyFilesFromCalcLoc(
                calc_loc=fw_outcar_to_born_name,
                filenames=['BORN']
            ),
            CopyFilesFromCalcLoc(
                calc_loc=fw_make_mask_name,
                filenames=['mask.yaml']
            ),
            # TODO: Replace the following task with one that writes to JSON
            PyTask(func='lblcrn.experiments.ir.wrapper.ir_sim_wrapper'),
            # TODO: Push it to the DB with MSONToDB
        ],
        name='IR Spectrum Calculation' + tracking_id,
        parents=[fw_outcar_to_born, fw_make_mask]
    )

    return [
        fw_vasp_relax,
        fw_vasp_born,
        fw_vasp_force_consts,
        fw_outcar_to_born,
        fw_make_mask,
        fw_ir_spectrum,
    ]