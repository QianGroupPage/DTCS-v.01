"""A built-in help function for the module lblcrn.

It's private so that from lblcrn import * does't override Python's built-in
help function.

Exposes:
    lblcrn_help: A help function, also accessible through lblcrn.help().
"""

import pathlib
import os
import webbrowser

import lblcrn


DOCS_PATH = lblcrn._DATA_FILES_PATH  + '/docs/index.html'
DOCS_URL = pathlib.Path(os.path.abspath(DOCS_PATH)).as_uri()


def lblcrn_help():
    """Helps you."""
    print('What do you want to be in the help function?')  # TODO


def lblcrn_docs():
    """Opens the docs in your web browser."""
    # TODO: Remove this later, we just don't have a docs host right now.
    webbrowser.open(DOCS_URL)
