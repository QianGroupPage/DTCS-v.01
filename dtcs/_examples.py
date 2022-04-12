"""Utility for packaged examples.

Exports:
    list_all()
    print_all()
    load(name: str)

Usage:
    list all of the example with print_all()::

        print_all()
        > Example 'predator-prey' (type RxnSystem)

    Then load one of them::

        rsys = load('predator-prey')
"""

import glob
import json
import os

import pkg_resources
import monty.json

import dtcs
from dtcs import _logger


EXAMPLE_OBJS_DIR = dtcs._DATA_FILES_PATH + 'examples/objects/'
EXAMPLE_OBJS_TYPE = '.json'


def list_all():
    """Return a list of all examples which can be loaded.

    `[load(example) for example in list_all()]` should work.
    """
    pattern = EXAMPLE_OBJS_DIR + '*' + EXAMPLE_OBJS_TYPE
    files = glob.glob(pattern)
    names = []
    for file in files:
        names.append(file.split(os.sep)[-1].replace(EXAMPLE_OBJS_TYPE, ''))
    return names


def print_all():
    """Print all of the examples which can be loaded."""
    for name in list_all():
        cls_ = json.loads(_read_example(name))['@class']
        print(f'Example \'{name}\' (type: {cls_})')


def load(name: str):
    """Loads the example given by name."""
    dct = _read_example(name)
    example = json.loads(dct, cls=monty.json.MontyDecoder)

    # Echo the example you're loading.
    _logger.echo(f'loaded: {example}')

    return example


def _read_example(name: str) -> str:
    """Loads the example given by name as a raw string.

    Args:
        name: The name of the example in the examples/objects directory.

    Returns:
        The str of the file examples/objects/{name}.json
    """
    # Get the path to the data file.
    ex_path = EXAMPLE_OBJS_DIR + name + EXAMPLE_OBJS_TYPE

    # We could actually just use open(), as ex_path is absolute,
    # but it might not be absolute later, so I'll leave it in.
    with pkg_resources.resource_stream(__name__, ex_path) as ex_file:
        example = ex_file.read()

    return example
