"""LBL Chemical Reaction Network Simulator

- Sets environment variable PYGAME_HIDE_SUPPORT_PROMPT to 'hide'
- Creates env variable LBLCRN_DO_ECHO
- Creates _DATA_FILES_PATH, path by which to access data files
- Has master __version__ copy
"""

import os
import sys
import warnings

from dtcs import config

# Hide the pygame support prompt before importing dtcs.
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = 'hide'

# A module-wide switch to (dis)allow printing through the _echo module.
# TODO: Every time you import the module it'll set this to false,
#  so we should edit it to check if the environment variable already exists.
#  so that it can be set before-hand.
os.environ['LBLCRN_DO_ECHO'] = 'false'

# The path for which to access data files.
# Needs to exist before importing dtcs._examples
_DATA_FILES_PATH = sys.prefix + '/dtcs/'

# Master version information. Modifying this should update everything else.
__version__ = '2023.11.19dev1'

# Enable depreciation warnings
warnings.filterwarnings(config.depr_warnings,
                        category=DeprecationWarning,
                        module='dtcs')

from dtcs._logger import dtcs_echo_on, dtcs_echo_off
import dtcs.resources as dtcs_resources
import dtcs.examples as dtcs_examples
from dtcs.help import *

# Make the following names fake so that `from dtcs import *` doesn't
# create mysterious (or conflicting) variables.
def __getattr__(name: str):
    if name == 'do_echo':
        if os.environ['LBLCRN_DO_ECHO'] == 'true':
            return True
        return False
    elif name == 'docs':
        return dtcs_docs
    elif name == 'echo_on':
        return dtcs_echo_on
    elif name == 'echo_off':
        return dtcs_echo_off
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


__all__ = [
    'dtcs_docs',
    'dtcs_echo_off',
    'dtcs_echo_on',
    'dtcs_docs',
    'dtcs_examples',
    'dtcs_help',
    'dtcs_resources',
]
