"""LBL Chemical Reaction Network Simulator

TODO: A docstring
"""

import os
import sys

# Hide the pygame support prompt before importing lblcrn.
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = 'hide'

# The path for which to access data files.
# Needs to exist before importing lblcrn._examples
_DATA_FILES_PATH = sys.prefix + '/lblcrn/'

# Master version information. Modifying this should update everything else.
__version__ = 'dev0.1.2.1'

# A module-wide switch to (dis)allow printing through the _echo module.
do_echo = False  # TODO: consider using os.environ?


from lblcrn.common import *
from lblcrn.crn_sym import *
from lblcrn.bulk_crn import *
from lblcrn.surface_crn import *
from lblcrn._echo import lblcrn_echo_on, lblcrn_echo_off
import lblcrn._examples as lblcrn_examples
from lblcrn._help import *

# Make the following names fake so that `from lblcrn import *` doesn't
# create mysterious (or conflicting) variables.
def __getattr__(name):
    if name == 'echo_on':
        return lblcrn_echo_on
    elif name == 'echo_off':
        return lblcrn_echo_off
    elif name == 'examples':
        return lblcrn_examples
    elif name == 'help':
        return lblcrn_help
    else:
        raise AttributeError(f'module \'{__name__}\' '
                             f'has no attribute \'{name}\'')
