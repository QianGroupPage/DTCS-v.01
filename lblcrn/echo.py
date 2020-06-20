"""Allow/disallow user prompts and messages.

Internally, this module is used for optional user-interaction, intended for
use in iPython notebooks.

Exposes:
    lblcrn_echo_on()
    lblcrn_echo_off()
"""

import lblcrn

def lblcrn_echo_on():
    """Turn on echo."""
    lblcrn.do_echo = True

def lblcrn_echo_off():
    """Turn off echo."""
    lblcrn.do_echo = False

def _input(*args, **kwargs):
    """Call input() if echo is on."""
    if lblcrn.do_echo:
        return input(*args, **kwargs)

def _input_yn(*args, **kwargs):
    """Call input() and interpret a yes/no answer."""
    ans = input(*args, **kwargs)
    if ans.lower() in ['y', 'ye', 'yes']:
        return True
    return False

def _print(*args, **kwargs):
    """Call print() if echo is on."""
    if lblcrn.do_echo:
        print(*args, **kwargs)
