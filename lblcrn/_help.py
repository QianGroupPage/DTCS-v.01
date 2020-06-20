"""A built-in help function for the module lblcrn.

It's private so that from lblcrn import * does't override Python's built-in
help function.

Exposes:
    lblcrn_help: A help function, also accessible through lblcrn.help().
"""

def lblcrn_help():
    """TODO
    """

    from lblcrn import echo
    print('TODO')  # TODO
