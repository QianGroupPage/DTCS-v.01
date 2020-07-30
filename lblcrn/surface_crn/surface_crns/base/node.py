class Node:
    '''
    Represents a spatial location. Has a chemical state represented by a
    string, a list of neighbors with weights on each connection, and a timestamp
    of the last update.

    self.neighbors is a list of 2-tuples, where the first element of each
    tuple is a neighboring node and the second element of each tuple is the
    weight of the edge connecting the two (defaults to 1, if no weight is
    given).
    '''
    def __init__(self, state = None, neighbors = None, timestamp = 0,\
                 position = None, group=None, node_id=None):
        if state == None:
            self.state = ""
        else:
            self.state = state
        if neighbors == None:
            self.neighbors = []
        elif isinstance(neighbors, list):
            def weighted_neighbor(neighbor):
                if isinstance(neighbor, tuple) and len(neighbor == 2):
                    return neighbor
                else:
                    return (neighbor, 1)
            self.neighbors = list(map(weighted_neighbor, neighbors))
        else:
            raise TypeError("Neighbors must be a list.")

        self.timestamp = timestamp
        if position == None:
            position = ()
        self.position = position

        # List of nodes in the same group, these should appear and disappear at the same time
        if group == None:
            self.group = []
        else:
            self.group = group

        self.node_id = node_id


    def is_empty(self):
        return self.state == ""

    def __str__(self):
        str_rep = "<State: " + self.state + " (updated at time " + \
            str(self.timestamp) + ")"
        if len(self.neighbors) > 0:
            str_rep += "; neighbor states: "
            str_rep += self.neighbors[0][0].state if self.neighbors[0] \
                        else "<node is NoneType>"
            for i in range(1, len(self.neighbors)):
                str_rep += ", "+self.neighbors[i][0].state if self.neighbors[i]\
                            else "<node is NoneType>"
        else:
            str_rep += "; No neighbors"

        if len(self.group) > 0:
            str_rep += "; neighbor states: "
            str_rep += self.group[0].state if self.group[0] \
                       else "<node is NoneType>"
            for i in range(1, len(self.group)):
                str_rep += ", "+self.group[i].state if self.group[i]\
                            else "<node is NoneType>"
        else:
            str_rep += "; No group"
        if self.position != ():
            str_rep += "; Position: " + str(self.position)

        str_rep += ">"
        return str_rep

    def __repr__(self):
        return str(self)
#end class Node
