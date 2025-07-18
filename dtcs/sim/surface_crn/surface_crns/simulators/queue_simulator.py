import math
import random
from queue import *

import numpy as np

from dtcs.sim.surface_crn.surface_crns.simulators.event import Event


class QueueSimulator:
    """
    Surface CRN simulator based on Gillespie-like next-reaction determination
    at each node. Upcoming reactions are stored in a priority queue, sorted
    by reaction time. Each time an event occurs, the next reaction time for
    each participating node is recalculated and added to the queue. Timestamps
    are used to ensure that a node does not react if it was changed between the
    time its reaction was issued and the time the reaction would occur.

    Uses unimolecular and bimolecular reactions only.
    """

    def __init__(self, surface=None, transition_rules=None, seed=None, group_selection_seed=None,
                 simulation_duration=100, rxns=None, sm=None):
        self.debugging = False
        if transition_rules is None:
            self.rule_set = []
        else:
            self.rule_set = transition_rules

        self.group_selection_rng = random.Random(group_selection_seed)
        # The primary pseudo-random number generator that decides the timestamps
        self.main_random_generator = random.Random(seed)

        self.simulation_duration = simulation_duration
        self.surface = surface
        self.sm = sm
        self.rxns = rxns
        self.add_groups()

        self.init_state = surface.get_global_state()
        # Build a mapping of states to the possible transitions they could
        # undergo.
        self.rules_by_state = dict()
        for rule in self.rule_set:
            for input_state in rule.inputs:
                if not input_state in self.rules_by_state:
                    self.rules_by_state[input_state] = []
                if not rule in self.rules_by_state[input_state]:
                    self.rules_by_state[input_state].append(rule)

        if self.rxns:
            self.rxns_by_rule_str = {}
            for rxn in self.rxns.surface_rxns:
                if rxn.is_reversible:
                    forward_rxn, backward_rxn = rxn.to_rxns()
                    self.rxns_by_rule_str[forward_rxn.surface_engine_str] = forward_rxn
                    self.rxns_by_rule_str[backward_rxn.surface_engine_str] = backward_rxn
                else:
                    self.rxns_by_rule_str[rxn.surface_engine_str] = rxn
            # self.surface_rxns_objects_by_rule = {rule: self.rxns_by_rule_str[str(rule)] for rule in self.rule_set}

        # TODO
        self.concentration_trajectory = None

        self.time = 0
        # self.surface.set_global_state(self.init_state)

        # TODO: Remove this hard-coded option, instead use logger
        if self.debugging:
            print("QueueSimulator initialized with global state:\n" + str(self.init_state))
        self.reset()
        if self.debugging:
            print(self.surface)

        self.initiate_concs()
        if self.surface:
            self.add_concs()
        # if rxns:
        #     self.initiate_marker_concs()

    def reset(self):
        """
        Clear any reactions in the queue and populate with available reactions.
        """
        self.event_queue = PriorityQueue()
        self.initialize_reactions()

    def initialize_reactions(self):
        """
        Populate the reaction queue with initial reactions.
        """
        for node in self.surface:
            node.timestamp = self.time
            self.add_next_reactions_with_node(node=node,
                                              first_reactant_only=True,
                                              exclusion_list=[])

    def done(self):
        """
        True iff there are no more reactions or the simulation has reached
        final time.
        """
        return self.event_queue.empty() or self.time >= self.simulation_duration

    def add_groups(self):
        """
        Update the surface with groups in the species manager.

        :param surface: a surface structure
        :param rsys: a rxn system
        """
        size_dict = self.sm.large_species_dict
        seen = set()
        for s in self.surface:
            if s not in seen and s.state in size_dict:
                group_size = size_dict[s.state]
                # Build the group
                free_neighbors = []
                for t in s.neighbors:
                    n = t[0]
                    if n not in seen and n.state in self.rxns.surface_names:
                        free_neighbors.append(n)
                    seen.add(n)
                seen.add(s)

                # TODO: if there are too few free_neighbors, do something else
                # Pick from free neighbors as part of the group
                group = self.group_selection_rng.sample(free_neighbors, group_size - 1) + [s]

                for n in group:
                    n.group = group
                    n.state = s.state

    def process_next_reaction(self):
        local_debugging = False
        """
        Process and return the next reaction in the queue:
        (1) Make sure the reaction is still valid (if not, try the next one
            instead).
        (2) Update the surface based on the reaction.
        (3) Determine the next reactions for each node involved in the reaction
            and add them to the event queue.
        """
        next_reaction = None
        while next_reaction == None:
            if self.event_queue.empty():
                self.time = self.simulation_duration
                return None

            # TODO: what's in "next_reaction"?
            next_reaction = self.event_queue.get()
            if next_reaction.time > self.simulation_duration:
                self.time = self.simulation_duration
                return None
            self.time = next_reaction.time
            participants = next_reaction.participants
            # TODO: add a group of species when appropriate
            outputs = next_reaction.rule.outputs
            if local_debugging:
                print("Processing event " + str(next_reaction.rule) +
                      " at time " + str(self.time) + ", position " +
                      str(participants[0].position) + " ")

            # # TODO: print out the surface
            # print("Processing event " + str(next_reaction.rule) +
            #       " at time " + str(self.time) + ", position " +
            #       str(participants[0].position) + " ")

            # If the first input was modified since the event was issued,
            # don't run it
            if participants[0].timestamp > next_reaction.time_issued:
                if local_debugging:
                    print("ignored -- first reactant changed since event " +
                          "issued.")
                next_reaction = None
                continue

            outgoing_group = []
            if len(participants) > 1:
                if local_debugging:
                    print("and " + str(participants[1].position) + " ")
                # If the second input was modified since the event was
                # issued, don't run it
                if participants[1].timestamp > next_reaction.time_issued:
                    if local_debugging:
                        print("ignored -- second reactant changed since " +
                              "event issued.")
                    next_reaction = None
                    continue
                # TODO: accomodations for 2-sized species
                self.update_node(participants[1], outputs[1], outgoing_group)
            # TODO: accomodations for 2-sized species
            self.update_node(participants[0], outputs[0], outgoing_group)

            for member_node in participants[0].group:
                if member_node is not participants[0]:
                    next_reaction.participants.append(member_node)
            if len(participants) > 1:
                for member_node in participants[1].group:
                    if member_node is not participants[1]:
                        next_reaction.participants.append(member_node)
            # Signify that the out-going group are also participants
            next_reaction.participants.extend(outgoing_group)

            if local_debugging:
                print("processed.")

            # Determine the next reactions performed by each participant
            # changed in this reaction.
            if local_debugging:
                print("Checking for new reactions with node:" + \
                      str(participants[0]))
            self.add_next_reactions_with_node(participants[0],
                                              first_reactant_only=False,
                                              exclusion_list=[])
            if len(participants) > 1:
                if local_debugging:
                    print("Checking for new reactions with node:" + \
                          str(participants[1]))
                # TODO: build up to the next reaction from here.
                self.add_next_reactions_with_node(
                    participants[1],
                    first_reactant_only=False,
                    exclusion_list=[participants[0]])
        if local_debugging:
            #     # TODO
            # if (participants[0].group and len(participants[0].group) > 1) or \
            #         (len(participants) > 1 and participants[1].group and len(participants[1].group) > 1):
            print("process_next_reaction() returning event " +
                  str(next_reaction))

        self.add_concs()
        # self.update_markers(next_reaction) TODO(Andrew)
        return next_reaction

    # end def process_next_reaction

    def update_node(self, node, new_state, outgoing_members):
        """
        This function updates node to the new_state.

        If node has a group, the group will be removed. If new_state has a group, randomly select node's free neighbors
        to form a new group.

        :param node: the node to be updated;
        :param new_state: the new state for the node;
        :param outgoing_members: a list to contain other members node's group, in case that the entire group should be
        removed.
        :return: None
        """
        # TODO
        # print("updating node")
        # print(f'{self.time}, from {node.state} at {node.position} to {new_state}')
        output_state = new_state
        default_state = self.sm.site_species_map(self.sm, output_state)
        if self.sm and output_state in self.sm.large_species_dict:
            free_neighbors = 0
            for neighbor_node, _ in node.neighbors:
                if neighbor_node.state == default_state:
                    free_neighbors += 1
            group_size = self.sm.large_species_dict[output_state]
            # No reaction would happen if not enough free neighbors
            if free_neighbors + 1 < group_size:
                # print("not enough free neighbors")
                return

        original_state = node.state
        default_state = self.sm.site_species_map(self.sm, original_state)

        if len(node.group) > 1:
            # print(f"The default state for {original_state} should be {default_state}.")

            for neighbor_node in node.group:
                if neighbor_node is not node:
                    # print("resetting neighbor node", neighbor_node.state)
                    neighbor_node.group = []
                    # Other members of the group should be set to default surface species.
                    neighbor_node.state = default_state
                    # This implies the node should not be used in this step anymore
                    neighbor_node.timestamp = self.time
                    outgoing_members.append(neighbor_node)
            node.group = []
            node.state = default_state

        if self.sm and output_state in self.sm.large_species_dict:
            new_group = [node]
            free_neighbors = []
            for neighbor_node, _ in node.neighbors:
                if neighbor_node.state == default_state:
                    free_neighbors.append(neighbor_node)
            group_size = self.sm.large_species_dict[output_state]
            # TODO: set the seed, to preserve simulation reproducibility.
            new_group += self.group_selection_rng.sample(free_neighbors, group_size - 1)

            # print("\n\n\n\n\n")
            # print("size of new group", len(new_group))
            for member_node in new_group:
                member_node.state = output_state
                member_node.timestamp = self.time
                member_node.group = new_group
        else:
            node.state = output_state
            node.timestamp = self.time

    # end def update_node

    # TODO:
    def add_concs(self):
        """
        Add the concentrations from the current timestamps to concentration_trajectory dataframe.
        """
        # TODO: sum things up here.
        # TODO: easy ways to average out the last 2s, including cutting out the "unqualified" seconds each
        # time before you do averaging.
        species_count = self.surface.species_count()

        for k, l in self.concs.items():
            if k not in species_count:
                l.append(0)
            else:
                l.append(species_count[k])

        for k, v in self.surface.species_count().items():
            if k not in self.concs:
                self.concs[k] = [0 for _ in self.times]
                self.concs[k].append(v)

        # self.concs.append(species_count)
        #
        self.times.append(self.time)
        #
        # current_row = pd.DataFrame(species_count, index=[self.time])
        #
        # if self.concentration_trajectory is None:
        #     self.concentration_trajectory = current_row
        # else:
        #     # TODO: this will cause serious speed issues. Fix them.
        #     self.concentration_trajectory = self.concentration_trajectory.append(current_row)
        # self.concentration_trajectory.fillna(0, inplace=True)

    def initiate_concs(self):
        self.concs = {}
        self.times = []

    def add_next_reactions_with_node(self, node, first_reactant_only=False,
                                     exclusion_list=None):
        """
        Determines whether or not the specified node can react according to any
        known rule and, if it can, returns a new event for the next occurrence
        of that reaction.

        If first_reactant_only = True, only checks to see if the species is the
        FIRST input species. For example, a node with state 'A' will react
        according to the rule "A + B -> C + D", but NOT according to the rule
        "B + A -> D + C". This mode is intended to make initialization of the
        surface easier.

        Nodes in exclusion_list are not considered eligible for reaction.
        """
        local_debugging = False
        if exclusion_list == None:
            exclusion_list = []
        if node.state not in self.rules_by_state:
            if local_debugging:
                print("No reactions possible with state " + node.state)
            return
        for rule in self.rules_by_state[node.state]:
            if local_debugging:
                print("Checking rule\n\t" + str(rule) + "\nagainst node\n\t" + \
                      str(node))
            if (first_reactant_only and node.state != rule.inputs[0]) or \
                    (not first_reactant_only and not node.state in rule.inputs):
                # Not eligible for reaction. Don't add.
                if local_debugging:
                    print("Rule rejected: input 1 doesn't match.")
                continue

            # If it's a unimolecular reaction, we can now add the reaction to
            # the event queue.
            if len(rule.inputs) == 1:
                # TODO: marker addition
                # TODO: this determines the time step till next reaction.
                time_to_reaction = np.log(1.0 / self.main_random_generator.random()) / rule.rate
                if math.isinf(time_to_reaction):
                    continue
                event_time = self.time + time_to_reaction
                new_event = Event(time=event_time,
                                  rule=rule,
                                  participants=[node],
                                  time_issued=self.time)
                self.event_queue.put(new_event)
                if local_debugging:
                    print("Event added: " + str(new_event))
                    print(str(new_event))

            elif len(rule.inputs) == 2:
                # If a bimolecular reaction, then we need to keep track of
                # which reactant this node is, so that we can check the other
                # one against each neighbor.
                if first_reactant_only:
                    node_index = 0
                else:
                    node_index = rule.inputs.index(node.state)
                if local_debugging:
                    print("First node has index " + str(node_index))

                for neighbor_node, weight in node.neighbors:
                    if local_debugging:
                        print("\tChecking neighbor node with state " + \
                              neighbor_node.state)
                    if neighbor_node.state != rule.inputs[1 - node_index]:
                        # Not eligible for reaction
                        if local_debugging:
                            print("Rule rejected: input 2 doesn't match.")
                        continue
                    if neighbor_node in exclusion_list:
                        if local_debugging:
                            print("Neighbor is on exclusion list. Moving on.")
                        continue

                    # If the two inputs are identical and the this node could
                    # be either reactant, then the reaction needs to be counted
                    # twice.
                    if not first_reactant_only and \
                            rule.inputs[0] == rule.inputs[1]:
                        num_reactions = 2
                    else:
                        num_reactions = 1
                    for x in range(num_reactions):
                        rate = rule.rate * weight
                        time_to_reaction = np.log(1.0 / self.main_random_generator.random()) / rate
                        if math.isinf(time_to_reaction):
                            continue
                        event_time = self.time + time_to_reaction
                        new_participants = [None, None]
                        new_participants[node_index] = node
                        new_participants[1 - node_index] = neighbor_node
                        # If counting the reaction twice, it needs to be flipped
                        # for the second reaction. Example: A + A -> B + C
                        if x == 2:
                            new_participants.reverse()
                        new_event = Event(time=event_time,
                                          rule=rule,
                                          participants=new_participants,
                                          time_issued=self.time)
                        self.event_queue.put(new_event)
                        if local_debugging:
                            print("Event added: " + str(new_event))
            else:
                raise Exception("Error in transition rule " + str(rule) + \
                                "\nOnly rules with one or two inputs allowed!")

    # end def add_next_reactions_with_node

    # Used only when RXN System is provided to the simulator.
    def update_markers(self, reaction):
        """
        This should run after add_concs is run.

        Update an existing marker_concs, already with initial values for all markers.
        """
        # In case only text-based manifest files are provided.
        if not self.rxns:
            return
        rxn = self.rxns_by_rule_str[reaction.rule.identifiable_string()]

        for conc in self.marker_concs.values():
            conc.append(conc[-1])
        for marker in rxn.markers.reactants:
            self.marker_concs[marker][-1] -= 1
        for marker in rxn.markers.products:
            self.marker_concs[marker][-1] += 1

    def initiate_marker_concs(self):
        # From a Marker object to its trajectory list.
        self.marker_concs = {marker: [marker.initial_count] for marker in self.rxns.species_manager.get_all_markers()}

    def drop_data(self):
        """
        Drop all trajectory data stored in this simulator
        :return:
        """
        for m, conc in self.marker_concs.items():
            m.initial_count = conc[-1] if conc else 0
        # self.initiate_marker_concs() TODO(Andrew)
        self.initiate_concs()
# end class QueueSimulator
