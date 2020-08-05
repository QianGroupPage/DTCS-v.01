"""VASP Fireworks for LBLCRN

TODO(Andrew)
"""

from typing import Optional

from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from fireworks import Firework
from pymatgen.io.vasp.sets import MPRelaxSet

from lblcrn.crn_sym.rxn_system import RxnSystem
from lblcrn.vasp.firetasks.parse_outputs import MsonToDb
from lblcrn.vasp.firetasks.crn import BulkCRNSim
from lblcrn.vasp.firetasks.xps import SimulateXPS


__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'


class CRNSimulateFW(Firework):  # TODO: Move this out of the vasp subpackage.
    """
    TODO
    """

    def __init__(
            self,
            reaction_system: RxnSystem,
            sim_type: str = 'bulk',
            sim_options: Optional[dict] = None,
            name: Optional[str] = None,
            db_file=None,
            wf_meta=None,
            **kwargs,
    ):
        """
        TODO
        """
        name = name or f'{sim_type} CRN simulation'
        sim_options = sim_options or {}

        if sim_type == 'bulk':
            crn_sim_task = BulkCRNSim(
                reaction_system=reaction_system.as_dict(),
                sim_options=sim_options,
            )
        else:
            raise ValueError(f'Invalid simulation type {sim_type}.')

        tasks = [
            crn_sim_task,
            PassCalcLocs(name=name),
            MsonToDb(
                in_fname='crn_time_series.json',
                db_name='crn_sims',
                db_file=db_file,
                wf_meta=wf_meta,
            ),
        ]

        super().__init__(
            tasks=tasks,
            name=name,
            **kwargs
        )


class XPSSimulateFW(Firework):  # TODO: Move this out of the vasp subpackage.
    """
    TODO
    """

    def __init__(
            self,
            name: str = None,
            db_file=None,
            wf_meta=None,
            **kwargs,
    ):
        """
        TODO
        """
        name = name or 'XPS simulation'

        xps_sim_task = SimulateXPS(
            # TODO
        )

        tasks = [
            xps_sim_task,
            #PassCalcLocs(name=name),
            #MsonToDb(
            #    in_fname='xps.json',
            #    db_name='xps',
            #    db_file=db_file,
            #    wf_meta=wf_meta,
            #),
        ]

        super().__init__(
            tasks=tasks,
            name=name,
            **kwargs
        )