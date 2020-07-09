import numpy as np
import pandas as pd
import random
import math
from queue import *
from lblcrn.surface_crn.surface_crns.simulators.event import Event

class QueueSimulator:
    '''
    Surface CRN simulator based on Gillespie-like next-reaction determination
    at each node. Upcoming reactions are stored in a priority queue, sorted
    by reaction time. Each time an event occurs, the next reaction time for
    each participating node is recalculated and added to the queue. Timestamps
    are used to ensure that a node does not react if it was changed between the
    time its reaction was issued and the time the reaction would occur.

    Uses unimolecular and bimolecular reactions only.
    '''
    def __init__(self, surface = None, transition_rules = None, seed = None,
                 simulation_duration = 100):
        self.debugging = False
        if transition_rules == None:
            self.rule_set = []
        else:
            self.rule_set = transition_rules

        random.seed(seed)
        self.simulation_duration = simulation_duration
        self.surface = surface
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

        # TODO
        self.concentration_trajectory = None

        self.time = 0
        self.surface.set_global_state(self.init_state)
        if self.debugging:
            print("QueueSimulator initialized with global state:\n" + str(self.init_state))
        self.reset()
        if self.debugging:
            print(self.surface)

    # TODO:
    def add_concs(self):
        """
        Add the concentrations from the current timestamps to concentration_trajectory dataframe.
        """
        species_count = dict(self.surface.species_count())

        if self.concentration_trajectory is None:
            self.concentration_trajectory = pd.DataFrame(species_count, index=[self.time])
        else:
            self.concentration_trajectory = self.concentration_trajectory.append(pd.DataFrame(species_count, index=[self.time]))

    def reset(self):
        '''
        Clear any reactions in the queue and populate with available reactions.
        '''
        self.event_queue = PriorityQueue()
        self.initialize_reactions()

    def initialize_reactions(self):
        '''
        Populate the reaction queue with initial reactions.
        '''
        for node in self.surface:
            node.timestamp = self.time
            self.add_next_reactions_with_node(node=node,
                                              first_reactant_only=True,
                                              exclusion_list = [])
    def done(self):
        '''
        True iff there are no more reactions or the simulation has reached
        final time.
        '''
        return self.event_queue.empty() or self.time >= self.simulation_duration

    def process_next_reaction(self):
        local_debugging = False
        '''
        Process and return the next reaction in the queue:
        (1) Make sure the reaction is still valid (if not, try the next one
            instead).
        (2) Update the surface based on the reaction.
        (3) Determine the next reactions for each node involved in the reaction
            and add them to the event queue.
        '''
        next_reaction = None
        while next_reaction == None:
            if self.event_queue.empty():
                self.time = self.simulation_duration
                return None

            next_reaction = self.event_queue.get()
            if next_reaction.time > self.simulation_duration:
                self.time = self.simulation_duration
                return None
            self.time     = next_reaction.time
            participants  = next_reaction.participants
            # TODO: add a group of species when appropriate
            outputs       = next_reaction.rule.outputs
            if local_debugging:
                print("Processing event " + str(next_reaction.rule) +
                      " at time " + str(self.time) + ", position " +
                      str(participants[0].position) + " ")

            # If the first input was modified since the event was issued,
            # don't run it
            if participants[0].timestamp > next_reaction.time_issued:
                if local_debugging:
                    print("ignored -- first reactant changed since event " +
                          "issued.")
                next_reaction = None
                continue

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
                # Change second reactant
                participants[1].state = outputs[1]
                participants[1].timestamp = self.time
            # Change first reactant
            participants[0].state = outputs[0]
            participants[0].timestamp = self.time

            # TODO: group edits:
            # Remove the outgoing group
            if len(participants[0].group) > 1:
                for n in participants[0].group:
                    n.group = []
                    # Other members of the group should be set to default surface species.
                    # n.state =
                    n.timestamp = self.time
                participants[0].group = []


            if local_debugging:
                print("processed.")

            # Determine the next reactions performed by each participant
            # changed in this reaction.
            if local_debugging:
                print("Checking for new reactions with node:" + \
                      str(participants[0]))
            self.add_next_reactions_with_node(participants[0],
                                              first_reactant_only = False,
                                              exclusion_list = [])
            if len(participants) > 1:
                if local_debugging:
                    print("Checking for new reactions with node:" + \
                          str(participants[1]))
                self.add_next_reactions_with_node(
                                            participants[1],
                                            first_reactant_only = False,
                                            exclusion_list = [participants[0]])
        if local_debugging:
            print("process_next_reaction() returning event " +
                  str(next_reaction))

        self.add_concs()
        return next_reaction
    #end def process_next_reaction

    def add_next_reactions_with_node(self, node, first_reactant_only = False,
                                        exclusion_list = None):
        '''
        Determines whether or not the specified node can react according to any
        known rule and, if it can, returns a new event for the next occurrence
        of that reaction.

        If first_reactant_only = True, only checks to see if the species is the
        FIRST input species. For example, a node with state 'A' will react
        according to the rule "A + B -> C + D", but NOT according to the rule
        "B + A -> D + C". This mode is intended to make initialization of the
        surface easier.

        Nodes in exclusion_list are not considered eligible for reaction.
        '''
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
            if (first_reactant_only     and node.state != rule.inputs[0]) or \
               (not first_reactant_only and not node.state in rule.inputs):
                    # Not eligible for reaction. Don't add.
                    if local_debugging:
                        print("Rule rejected: input 1 doesn't match.")
                    continue

            # If it's a unimolecular reaction, we can now add the reaction to
            # the event queue.
            if len(rule.inputs) == 1:
                time_to_reaction = np.log(1.0 / random.random()) / rule.rate
                if math.isinf(time_to_reaction):
                    continue
                event_time       = self.time + time_to_reaction
                new_event = Event(time = event_time,
                                  rule = rule,
                                  participants = [node],
                                  time_issued = self.time)
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
                    if neighbor_node.state != rule.inputs[1-node_index]:
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
                        time_to_reaction = np.log(1.0/random.random())/rate
                        if math.isinf(time_to_reaction):
                            continue
                        event_time       = self.time + time_to_reaction
                        new_participants = [None, None]
                        new_participants[node_index] = node
                        new_participants[1-node_index] = neighbor_node
                        # If counting the reaction twice, it needs to be flipped
                        # for the second reaction. Example: A + A -> B + C
                        if x == 2:
                            new_participants.reverse()
                        new_event = Event(time = event_time,
                                          rule = rule,
                                          participants = new_participants,
                                          time_issued = self.time)
                        self.event_queue.put(new_event)
                        if local_debugging:
                            print("Event added: " + str(new_event))
            else:
                raise Exception("Error in transition rule " + str(rule) + \
                                "\nOnly rules with one or two inputs allowed!")
    #end def add_next_reactions_with_node
# end class QueueSimulator