"""VASP Fireworks for LBLCRN

TODO(Andrew)
"""

from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.utils.utils import env_chk, get_logger
from fireworks import Firework

from lblcrn.vasp.firetasks.xps import SimulateXPS
from lblcrn.vasp.firetasks.parse_outputs import MsonToDb

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'

_logger = get_logger(__name__)


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
            # TODO: Give the user the ability to use these?
            #  it would have to be written as like a forced_input = ..., as
            #  the use of such a parameter will invalidate anything pushed
            #  into the spec (unless I add the functionality to search _tasks).
        )

        tasks = [
            xps_sim_task,
            PassCalcLocs(name=name),
            MsonToDb(
                in_fname='xps.json',
                db_name='xps',
                db_file=db_file,
                wf_meta=wf_meta,
            ),
        ]

        super().__init__(
            tasks=tasks,
            name=name,
            **kwargs
        )