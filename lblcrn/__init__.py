"""LBL Chemical Reaction Network Simulator

TODO: A docstring
"""

from lblcrn.common import *
from lblcrn.crn_sym import *
from lblcrn.bulk_crn import *
from lblcrn.echo import *
from lblcrn._help import *

__version__ = 'dev0.1.2.1'
do_echo = False

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
