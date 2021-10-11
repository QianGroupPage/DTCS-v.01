from lblcrn.sim.surface_crn.surface_crns.readers.statements import *

'''
For external use
'''


def read_totalistic_rules(filename):
    '''
    Read the update rule for a totalistic cellular automata from a totalistic
    rule file. The rule is written as a series of update options of the form:

    neighbor_total, current_state -> new_state

    where any spaces are optional.

    The file may optionally be terminated by a line of the string
    "!END_TOTALISTIC_RULE"
    '''
    with open(filename, 'rU') as rule_file:
        return parse_totalistic_rule_stream(rule_file)


def parse_totalistic_rule_stream(rules_stream):
    '''
    Parse a stream of strings containing update options, as from a totalistic
    rule file.
    See documentation for read_totalistic_rules for a description of the
    totalistic update option format.
    '''
    update_options = dict()
    for line in rules_stream:
        if line.startswith(section_ends['totalistic_rule']):
            break
        if not (line.startswith("#") or
                line.startswith("%") or
                line.strip() == ""):
            parse_option(update_options, line)

    def update_rule(neighbor_states, self_state):
        neighbor_sum = sum(map(lambda s: int(s), neighbor_states))
        return update_options[(neighbor_sum, self_state)]

    return update_rule


'''
Internal use only
'''


def parse_option(update_options, line):
    # Identify neighbor count
    line_parts = list(map(lambda s: s.strip(), line.split(',')))
    try:
        neighbor_count = int(line_parts[0])
    except ValueError:
        raise Exception('Illegal neighbor count ' + line_parts[0] +
                        ' in update option "' + line + '": neighbor count ' +
                        'must be an integer.')

    # Identify self-state
    line = line_parts[1]
    line_parts = list(map(lambda s: s.strip(), line.split('->')))
    if len(line_parts) != 2:
        raise Exception('Invalid update option "' + line + '": Input and ' +
                        'output states must be separated by a single "->".')
    self_state = line_parts[0]

    # Identify the new state
    new_state = line_parts[1]

    # Add the new rule
    update_options[(neighbor_count, self_state)] = new_state
