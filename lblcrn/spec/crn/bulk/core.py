"""TODO"""


from lblcrn.spec.crn.core import CRNSpec
from lblcrn.spec.species import SpeciesCollection
from lblcrn.spec.crn.bulk.rxn_system import BulkRxnSystem


class BulkCRNSpec(CRNSpec):
    """TODO
    """

    _default = {
        'sim_type': 'bulk',
        'time': 10,
        'max_step': 0.01,
    }