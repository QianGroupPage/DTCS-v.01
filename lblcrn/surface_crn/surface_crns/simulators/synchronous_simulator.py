import random
from lblcrn.surface_crn.surface_crns.base.transition_rule import TransitionRule
from lblcrn.surface_crn.surface_crns.simulators.event import Event

class SynchronousSimulator:
    '''
    Simulator for synchronous cellular automata. Updates according to
    update_rule, which is a function taking parameters "neighbor_states" and
    "current_state" and returning the next state for that node, based on the
    neighbors' states and the node's own state.
    '''
    def __init__(self, surface = None, update_rule = None, seed = None,
                 simulation_duration = 100):
        self.debugging = False
        if update_rule == None:
            self.update_rule = []
        else:
            self.update_rule = update_rule

        random.seed(seed)
        self.simulation_duration = simulation_duration
        self.surface = surface
        self.init_state = surface.get_global_state()

        self.initialize()

    def initialize(self):
        '''
        Start the simulation from the initial condition.
        '''
        self.time = 0
        self.surface.set_global_state(self.init_state)

    def done(self):
        '''
        True iff the simulation has reached final time.
        '''
        return self.time >= self.simulation_duration

    def process_next_reaction(self):
        local_debugging = False
        '''
        Update the surface one clock tick according to the rule given in
        update_rule
        '''
        changed_nodes = []
        for node in self.surface:
            neighbor_states = list(map(lambda tup:tup[0].state, node.neighbors))
            new_state = self.update_rule(neighbor_states, node.state)
            if node.state != new_state:
                changed_nodes.append(node)
                node.new_state = new_state

        # Have to make a second pass, because we can't make any changes until
        # all nodes have been checked.
        for node in changed_nodes:
            node.state = node.new_state

        # Return a bogus Event object to tell the controller which nodes to
        # update
        rule = TransitionRule(inputs = [0]*len(changed_nodes),
                              outputs = [1]*len(changed_nodes))
        new_event = Event(time = self.time,
                          rule = rule,
                          participants = changed_nodes,
                          time_issued = self.time)
        self.time += 1
        return new_event
    #end def process_next_reaction
# end class SynchronousSimulator
