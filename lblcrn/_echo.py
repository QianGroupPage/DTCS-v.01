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

# --- For Internal Use -------------------------------------------------------

def prompt(msg, *args, **kwargs):
    """Call input() if echo is on."""
    if lblcrn.do_echo:
        return input(msg, *args, **kwargs)

def prompt_yn(msg):
    """Call input() and interpret a yes/no answer.

    This will add ' (y/n): ' to your message, so just ask it a question like
    'The experimental data has a gaussian at 123 eV, mark as contaminant?'
    """
    if lblcrn.do_echo:
        ans = input(msg + ' (y/n): ')
        if ans.lower() in ['y', 'ye', 'yes']:
            return True
    return False

def echo(*args, **kwargs):
    """Call print() if echo is on."""
    if lblcrn.do_echo:
        print(*args, **kwargs)
