"""VASP Fireworks for LBLCRN

TODO(Andrew)
"""

from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from fireworks import Firework
from pymatgen.io.vasp.sets import MPRelaxSet

from lblcrn.vasp.firetasks.parse_outputs import CoreEigenToDb


__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'


class CRNSimulateFW(Firework):
    """
    TODO
    """

    def __init__(
            self,
            name=None,
            sim_kind: str = 'bulk',
            **kwargs,
    ):
        """
        TODO
        """

        tasks = []

        super().__init__(
            tasks=tasks,
            name=name,
            **kwargs
        )
