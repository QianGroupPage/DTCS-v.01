"""TODO"""

from typing import Dict, List, Mapping, Set

import sympy as sym

from lblcrn.spec.crn.crn_abc import CRNSpecABC
from lblcrn.spec.crn.bulk.rxn_system import BulkRxnSystem
from lblcrn.spec.species import SpeciesManager


class BulkCRNSpec(CRNSpecABC):
    """TODO
    """

    _rxn_sys_cls = BulkRxnSystem

    def __init__(self,
                 *components,
                 time: int = 10,
                 max_step: float = 0.01,
                 rsys: BulkRxnSystem = None,
                 species: SpeciesManager = None,
                 sim_type: str = None,
                 **kwargs):
        super().__init__(*components,
                         rsys=rsys,
                         species=species,
                         sim_type='bulk')
        self.time = time
        self.max_step = 0.1
