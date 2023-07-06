"""A built-in help function for the module dtcs.

It's private so that from dtcs import * does't override Python's built-in
help function.

Exposes:
    dtcs_help: A help function, also accessible through dtcs.help().
"""

import pathlib
import os
import webbrowser

import dtcs


DOCS_PATH = dtcs._DATA_FILES_PATH  + '/docs/index.html'
DOCS_URL = pathlib.Path(os.path.abspath(DOCS_PATH)).as_uri()


def dtcs_help():
    """Helps you."""
    print('What do you want to be in the help function?')  # TODO


def dtcs_docs():
    """Opens the docs in your web browser."""
    # TODO: Remove this later, we just don't have a docs host right now.
    webbrowser.open(DOCS_URL)
