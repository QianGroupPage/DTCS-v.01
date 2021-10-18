#!/usr/bin/env python
"""
lblcrn.__main__.py
------------------

An entry-point for the package for the purpose of help and basic tests: runs when
lblcrn is run as a module with `python -m lblcrn`.
"""

def cmd_help():
    """Show the help string."""

    help = """
    Basic CLI for lblcrn; for more information, refer to the docs and tutorial 
    at <https://github.com/AndrewBogdan/lbl-crn> or run lblcrn.lblcrn_help().
    
    Usage: `python -m lblcrn [options]`

    Options:
        -h, --help: Help
        -t, --test: Test to verify successful installation
        -V, --version: Show version number
    """
    print(help)


def cmd_version():
    """Gives the version of the package."""

    import lblcrn
    print(f'lblcrn version {lblcrn.__version__}')


if __name__ == '__main__':
    import sys

    opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
    args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]

    if '-h' in opts or '--help' in opts:
        cmd_help()
    elif '-V' in opts or '--version' in opts:
        cmd_version()
    else:
        cmd_help()
