"""LBL Chemical Reaction Network Simulator

- Sets environment variable PYGAME_HIDE_SUPPORT_PROMPT to 'hide'
- Creates env variable LBLCRN_DO_ECHO
- Creates _DATA_FILES_PATH, path by which to access data files
- Has master __version__ copy
"""

import os
import sys

# Hide the pygame support prompt before importing lblcrn.
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = 'hide'

# A module-wide switch to (dis)allow printing through the _echo module.
# TODO: Every time you import the module it'll set this to false,
#  so we should edit it to check if the environment variable already exists.
#  so that it can be set before-hand.
os.environ['LBLCRN_DO_ECHO'] = 'false'

# The path for which to access data files.
# Needs to exist before importing lblcrn._examples
_DATA_FILES_PATH = sys.prefix + '/lblcrn/'

# Master version information. Modifying this should update everything else.
__version__ = 'dev0.1.2.3'

# TODO(Andrew): I removed the import *s because importing everything was creating
#  issues with conflicts and syntax-error-y in-development code. Add it back
#  when you merge master.
from lblcrn._echo import lblcrn_echo_on, lblcrn_echo_off
import lblcrn._resources as lblcrn_resources
import lblcrn._examples as lblcrn_examples
from lblcrn._help import *

# Make the following names fake so that `from lblcrn import *` doesn't
# create mysterious (or conflicting) variables.
def __getattr__(name: str):
    if name == 'do_echo':
        if os.environ['LBLCRN_DO_ECHO'] == 'true':
            return True
        return False
    elif name == 'docs':
        return lblcrn_docs
    elif name == 'echo_on':
        return lblcrn_echo_on
    elif name == 'echo_off':
        return lblcrn_echo_off
    elif name == 'examples':
        return lblcrn_examples
    elif name == 'help':
        return lblcrn_help
    elif name == 'resources':
        return lblcrn_resources
    else:
        raise AttributeError(f'module \'{__name__}\' '
                             f'has no attribute \'{name}\'')

def __setattr__(name: str, value):
    if name == 'do_echo':
        if value:
            os.environ['LBLCRN_DO_ECHO'] = 'true'
        else:
            os.environ['LBLCRN_DO_ECHO'] = 'false'
    else:
        super().__setattr__(name, value)