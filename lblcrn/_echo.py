"""Allow/disallow user prompts and messages.

Internally, this module is used for optional user-interaction, intended for
use in iPython notebooks.

Exposes:
    lblcrn_echo_on()
    lblcrn_echo_off()

TODO: There _has_ to be a better way to do this, perhaps with the std
 logger and warnings module. I'll fix this eventually, hopefully - Andrew
"""

import lblcrn

@property
def do_echo():
    """Return lblcrn.do_echo"""
    return lblcrn.do_echo

def lblcrn_echo_on():
    """Turn on echo."""
    return EchoContext(True)


def lblcrn_echo_off():
    """Turn off echo."""
    return EchoContext(False)


# --- For Internal Use -------------------------------------------------------

class EchoContext:
    """
    TODO: I have no idea if this is the Pythonic way to do this, but I need it
     for debugging, so... yeah. - Andrew (andrewbogdan@lbl.gov)
    """

    def __init__(self, do_echo):
        self.do_echo = do_echo
        self.old_echo = lblcrn.do_echo

        lblcrn.do_echo = self.do_echo

    def __enter__(self):
        """We turn in echo in __init__, not in __enter__, as then you don't
        have to use a with statement to turn on echo."""
        pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        lblcrn.do_echo = self.old_echo


def prompt(msg, *args, **kwargs):
    """Call input() if echo is on."""
    if lblcrn.do_echo:
        return input(msg, *args, **kwargs)


def prompt_yn(msg, default=False):
    """Call input() and interpret a yes/no answer.

    This will add ' (y/n): ' to your message, so just ask it a question like
    'The experimental data has a gaussian at 123 eV, mark as contaminant?'
    """
    if lblcrn.do_echo:
        ans = input(msg + ' (y/n): ')
        if ans.lower() in ['y', 'ye', 'yes']:
            return True
        elif ans.lower() in ['n', 'no']:
            return False
    return default


def echo(*args, **kwargs):
    """Call print() if echo is on."""
    if lblcrn.do_echo:
        print(*args, **kwargs)
