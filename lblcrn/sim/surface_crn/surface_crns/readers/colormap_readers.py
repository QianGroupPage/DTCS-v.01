from lblcrn.sim.surface_crn.surface_crns.constants import *
from lblcrn.sim.surface_crn.surface_crns.ordered_dict import OrderedDict
from lblcrn.sim.surface_crn.surface_crns.readers.statements import *


def read_colormap(filename):
    """
    Read a file defining colors for state display. Definitions are written with
    one color on each line, in the following format:

    {state_class} state1, state2, state3: (r, g, b)

    where state1, state2, etc are names of a chemical states (as defined in the
    docs for transition rules) and r, g, b are values for red, green and blue
    between 0 and 255 (inclusive). Optionally, a group of states with the same
    color can be grouped into a class with name {state_class}.

    Lines beginning in '#' or '%' are comments.

    The colormap file may optionally be ended by a line starting with the string
    "!END_COLORMAP"
     """
    with open(filename, 'rU') as colormap_file:
        return parse_colormap_stream(colormap_file)


def parse_colormap_stream(colormap_stream):
    """
    Read a colormap from a stream of strings, as might be contained in a
    colormap definition file.
    See documentation for read_colormap for a description of the colormap file
    format.
    """
    colormap = OrderedDict()
    colormap[COLOR_CLASSES] = OrderedDict()
    for line in colormap_stream:
        if line.startswith(section_ends['colormap']):
            break
        if line.startswith('#') or line.startswith('%') or line.strip() == "":
            continue

        # Check for a class label
        class_label = None
        open_bracket_loc = line.find("{")
        if open_bracket_loc >= 0:
            close_bracket_loc = line.find("}")
            if close_bracket_loc < 0 or close_bracket_loc < open_bracket_loc:
                raise Exception('Invalid class declaration in colormap ' +
                                'definition: Missing or misplaced closing ' +
                                'bracket.')
            class_label = line[open_bracket_loc + 1:close_bracket_loc]
            line = line[:open_bracket_loc] + line[close_bracket_loc + 1:]

        # Check that there's a list of states and a color
        line_parts = line.split(":")
        if len(line_parts) != 2:
            raise Exception('Invalid color mapping definition "' + line +
                            '": State and color definition ' +
                            ' must be separated by a ":".')

        states = line_parts[0].split(',')
        states = list(map(lambda s: s.strip(), states))
        if class_label:
            if class_label in colormap[COLOR_CLASSES]:
                raise Exception('Two state classes cannot have the same name.')
            colormap[COLOR_CLASSES][class_label] = states
        color = line_parts[1].strip().strip('()')
        color_bits = color.split(',')
        color_bits = list(map(lambda s: int(s), color_bits))
        for state in states:
            if state == COLOR_CLASSES:
                raise Exception("Sorry, state string " +
                                COLOR_CLASSES +
                                " is reserved. Pick another state name.")
            colormap[state] = color_bits
    return colormap
