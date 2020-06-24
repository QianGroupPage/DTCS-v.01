"""LBL Chemical Reaction Network Simulator

TODO: A docstring
"""

from lblcrn.common import *
from lblcrn.crn_sym import *
from lblcrn.bulk_crn import *
from lblcrn.surface_crn import *
from lblcrn._echo import lblcrn_echo_on, lblcrn_echo_off
from lblcrn._help import *

__version__ = 'dev0.1.2.1'

# A module-wide switch to (dis)allow printing through the _echo module.
do_echo = False

# Make the following names fake so that `from lblcrn import *` doesn't
# create mysterious (or conflicting) variables.
def __getattr__(name):
    if name == 'help':
        return lblcrn_help
    elif name == 'echo_on':
        return lblcrn_echo_on
    elif name == 'echo_off':
        return lblcrn_echo_off
    else:
        raise AttributeError(f'module \'{__name__}\' '
                             f'has no attribute \'{name}\'')
