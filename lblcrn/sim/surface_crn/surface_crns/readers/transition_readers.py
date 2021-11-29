import re

from lblcrn.sim.surface_crn.surface_crns.base.transition_rule import TransitionRule
from lblcrn.sim.surface_crn.surface_crns.readers.statements import section_ends


"""
For external use
"""


def read_transition_rules(filename):
    """
    Read transition rules from a transition rule file. Transition rules are
    written with one rule on each line, in the following format:

    name1 + name2 -> name3 + name4 (<rate>)
    or
    name1 + name2 -> name3 + name4 {<rate>}

    where name1, etc are any alphanumeric labels for chemical states, and <rate>
    is a number. The rate may appear anywhere in the line, as long as it is
    contained within parentheses or curly brackets.

    Lines beginning with '#' or '%' are treated as comments.

    The file may optionally be terminated by a line of the string
    "!END_TRANSITION_RULES"
    """
    with open(filename, 'rU') as rule_file:
        return parse_transition_rule_stream(rule_file)


def parse_transition_rule_stream(rules_stream):
    """
    Parse a stream of strings containing  transition rules, as from a transition
    rule file.
    See documentation for read_transition_rules for a description of the
    transition rule format.
    """
    transition_rules = []
    for line in rules_stream:
        if line.startswith(section_ends['transition_rules']):
            break
        if not (line.startswith("#") or
                line.startswith("%") or
                line.strip() == ""):
            transition_rules.append(parse_rule(line))
    return transition_rules


"""
Internal use only
"""


def parse_rule(line):
    # Extract the rate constant
    open_paren_loc = line.find("(")
    if open_paren_loc < 0:
        open_paren_loc = line.find("{")
        close_paren_loc = line.find("}")
    else:
        close_paren_loc = line.find(")")
    if open_paren_loc < 0 or close_paren_loc < 0 or \
            close_paren_loc < open_paren_loc:
        raise Exception('Invalid transition rule "' + line +
                        '": No properly formatted rate constant.')
    rate_string = line[open_paren_loc + 1:close_paren_loc]
    try:
        slash_loc = rate_string.find('/')
        if slash_loc < 0:
            rate = float(rate_string)
        else:
            rate_numerator = float(rate_string[0:slash_loc])
            rate_denominator = float(rate_string[slash_loc + 1:])
            rate = rate_numerator / rate_denominator
    except ValueError:
        raise Exception('Invalid transition rule "' + line +
                        '": Improperly formatted rate constant ' + rate_string)
    line = line[:open_paren_loc] + line[close_paren_loc + 1:]

    # Split into input and output
    line_parts = line.split('->')
    if len(line_parts) != 2:
        raise Exception('Invalid transition rule "' + line + \
                        '": Input and output species must be ' + \
                        'separated by "->"')

    # Parse input species
    input_species = []
    input_strings = line_parts[0].split('+')
    for input_string in input_strings:
        try:
            input_string = parse_species(input_string)
            input_species.append(input_string)
        except SyntaxError:
            raise Exception('Invalid transition rule "' + line + \
                            '":Species ' + input_string + \
                            ' must be alphanumeric, but is not.')

    # Parse output species
    output_species = []
    output_strings = line_parts[1].split('+')
    for output_string in output_strings:
        try:
            output_string = parse_species(output_string)
            output_species.append(output_string)
        except SyntaxError:
            raise Exception('Invalid transition rule "' + line + \
                            '": Species ' + output_string + \
                            ' must be alphanumeric, but is not.')

    # Generate transition rule
    return TransitionRule(input_species, output_species, rate)


def parse_species(species_string):
    species_string = species_string.strip()
    # TODO: verif
    if not re.match('^[,_a-zA-Z0-9\*]+$', species_string):
        raise SyntaxError()
    return species_string
