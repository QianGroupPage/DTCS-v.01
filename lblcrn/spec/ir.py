"""TODO
"""

from lblcrn.spec import Spec
from lblcrn.spec.species import Species

__author__ = 'Andrew Bogdan'
__email__ = 'andrewbogdan@lbl.gov'


class IRSpecies(Species):
    """TODO"""

    _schema = [
        'vis_relax',
        'vis_born',
        'vis_fconsts',
    ]
