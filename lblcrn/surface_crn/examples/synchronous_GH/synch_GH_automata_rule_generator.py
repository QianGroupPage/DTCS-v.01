import itertools

direction_pairs = {'R':'L',
                   'U':'D',
                   'L':'R',
                   'D':'U'}
direction_transitions = {'U':'R',
                         'L':'U',
                         'D':'L',
                         'R':'D'}
colormap = {'Q':'(240, 240, 240)',
            'A':'(44,127,184)',
            'R':'(127,205,187)',
            'edge':'(160,160,160)'}
all_states = ['Q', 'A', 'R']
def state_pairs(): 
    return itertools.product(all_states, all_states)

def main():
    '''
    node state definition:
                <self_state>_<direction>_<position>_<states_seen>
    edge state definition:
                Edge_<direction> (Edges are non-activating)
    '''
    
    colormap_file = open("GH_spinning_arrow_colormap.txt",'w')

    # Write colormap rule for edge nodes separately.
    colormap_file.write('{Edge} ')
    first_color = True
    for direction in direction_pairs.keys():
        if first_color:
            first_color = False
        else:
            colormap_file.write(', ')
        colormap_file.write('Edge_' + direction)
    colormap_file.write(f": {colormap['edge']}\n")

    # Colormap rules for non-edges
    for state in all_states:
        colormap_file.write("{" + state + "} ")
        first_state = True
        for direction in direction_pairs.keys():
            for position in range(1,10):
                if first_state:
                    first_state = False
                else:
                    colormap_file.write(", ")
                colormap_file.write(f"{state}_{direction}_{position}_None")
                for s1 in all_states:
                    colormap_file.write(f", {state}_{direction}_{position}_{s1}")
                    for s2 in all_states:
                        colormap_file.write(f", {state}_{direction}_{position}_"
                                            f"{s1+s2}")
                        for s3 in all_states:
                            colormap_file.write(f", {state}_{direction}_"
                                                f"{position}_{s1+s2+s3}")
                            for s4 in all_states:
                                colormap_file.write(f", {state}_{direction}_"
                                                    f"{position}_{s1+s2+s3+s4}")
        colormap_file.write(f": {colormap[state]}\n")

    colormap_file.close()

    transition_rule_file = \
            open("GH_spinning_arrow_rules.txt", 'w')

    # Spinning transition rules, no edges
    for self_state, neighbor_state in state_pairs():
        for position in range(1,10):
            rule = spinning_rule("U", 
                                 self_state, 
                                 position, 
                                 "None", 
                                 neighbor_state, 
                                 "None")
            transition_rule_file.write(rule)
            for s1, ns1 in state_pairs():
                rule = spinning_rule("R", 
                                     self_state, 
                                     position, 
                                     s1, 
                                     neighbor_state, 
                                     ns1)
                transition_rule_file.write(rule)
                for s2, ns2 in state_pairs():
                    rule = spinning_rule("D", 
                                         self_state, 
                                         position, 
                                         s1+s2, 
                                         neighbor_state, 
                                         ns1+ns2)
                    transition_rule_file.write(rule)
                    for s3, ns3 in state_pairs():
                        rule = spinning_rule("R", 
                                             self_state, 
                                             position, 
                                             s1+s2+s3, 
                                             neighbor_state, 
                                             ns1+ns2+ns3)
                        transition_rule_file.write(rule)

    # Spinning transition rules against edges
    # Get the reset rules in here too, in the innermost loop.
    for self_state in all_states:
        for position in range(1,10):
            rule = edge_spinning_rule("U", 
                                      self_state, 
                                      position, 
                                      "None")
            transition_rule_file.write(rule)
            for s1 in all_states:
                rule = edge_spinning_rule("R", 
                                          self_state, 
                                          position, 
                                          s1)
                transition_rule_file.write(rule)
                for s2 in all_states:
                    rule = edge_spinning_rule("D", 
                                              self_state, 
                                              position, 
                                              s1+s2)
                    transition_rule_file.write(rule)
                    for s3 in all_states:
                        rule = edge_spinning_rule("R", 
                                                  self_state, 
                                                  position, 
                                                  s1+s2+s3)
                        transition_rule_file.write(rule)
                        # Get the reset rules here.
                        for s4 in all_states:
                            for direction in ["U", "D"]:
                                rule = reset_rule(direction,
                                                  self_state,
                                                  position,
                                                  s1+s2+s3+s4)
                                transition_rule_file.write(rule)

    # Close out transition rule file.
    transition_rule_file.close()


def spinning_rule(direction, self_state, position, states_seen, 
                           n_state, n_states_seen):
    n_direction = direction_pairs[direction]
    n_position = find_neighbor_position(direction, position)
    return (f"(1) {self_state}_{direction}_{position}_{states_seen} + "
           f"{n_state}_{n_direction}_{n_position}_{n_states_seen} -> "
           f"{self_state}_{direction_transitions[direction]}_{position}"
                f"_{('' if states_seen=='None' else states_seen)+n_state} + "
           f"{n_state}_{direction_transitions[n_direction]}_{n_position}"
                f"_{('' if n_states_seen=='None' else n_states_seen)+n_state}\n")

def edge_spinning_rule(direction, self_state, position, states_seen):
    edge_direction = direction_pairs[direction]
    return (f"(1) {self_state}_{direction}_{position}_{states_seen} + "
           f"Edge_{edge_direction} -> "
           f"{self_state}_{direction_transitions[direction]}_{position}"
                f"_{('' if states_seen=='None' else states_seen)}Q + Edge_{edge_direction}\n")

def reset_rule(direction, self_state, position, states_seen):
    if self_state == "A":
        new_state = "R"
    elif self_state == "R":
        new_state = "Q"
    elif self_state == "Q":
        if "A" in states_seen:
            new_state = "A"
        else:
            new_state = "Q"
    else:
        raise ValueError(f"Invalid emulated state {self_state}")
    return (f"(1) {self_state}_{direction}_{position}_{states_seen} -> "
            f"{new_state}_{direction}_{position}_None\n")


def find_neighbor_position(direction, position):
    pos_num = position - 1
    row_num = pos_num//3
    col_num = pos_num%3
    if direction == "U":
        row_num = (row_num-1)%3
    elif direction == "R":
        col_num = (col_num+1)%3
    elif direction == "D":
        row_num = (row_num+1)%3
    elif direction == "L":
        col_num = (col_num-1)%3
    pos_num = row_num*3 + col_num
    return pos_num + 1


if __name__ == "__main__":
    main()