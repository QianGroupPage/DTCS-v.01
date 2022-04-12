"""CRN Fireworks for LBLCRN

TODO(Andrew)
"""

from typing import Optional

from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.utils.utils import get_logger
from fireworks import Firework

from dtcs.spec.crn.crn_eni import RxnSystem
from dtcs.fire.firetasks.crn import BulkCRNSim
from dtcs.fire.firetasks.parse_outputs import MsonToDb

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'

_logger = get_logger(__name__)


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