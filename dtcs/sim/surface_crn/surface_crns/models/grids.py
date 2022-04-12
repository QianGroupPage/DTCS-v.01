from collections import Counter
from warnings import *

import numpy as np

from dtcs.sim.surface_crn.surface_crns.base.node import Node


class SquareGrid(object):
    """
    Representation of a CRN on a 2D finite rectangular mesh grid.

    Only allows reactions between directly-adjacent locations, and every
    edge has equal weight.

    Params:
        x_size, y_size: The number of cells in the x and y dimensions,
                        respectively.
        wrap: Iff true, the grid will wrap top to bottom and left to right.
                Default false.
    """

    def __init__(self, x_size, y_size, wrap=False):
        if not isinstance(x_size, int) or not isinstance(y_size, int):
            raise TypeError("SimGrid dimensions must be integers")
        self.x_size = x_size
        self.y_size = y_size
        self.wrap = wrap
        self.grid = np.empty([x_size, y_size], np.dtype(object))
        self.grid = np.array(self.grid)
        # Instantiate nodes
        self.populate_grid()

    def populate_grid(self):
        """
        Set/reset all nodes to new nodes with no state and populate neighbor
        nodes. Should only be used by initialization and set routines.
        """
        for x in range(self.x_size):
            for y in range(self.y_size):
                self.grid[x, y] = Node()
                self.grid[x, y].position = (x, y)
        # Populate node neighbor lists
        for x in range(self.x_size):
            for y in range(self.y_size):
                for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    nx = x + dx
                    ny = y + dy
                    if nx >= 0 and nx < self.x_size and \
                            ny >= 0 and ny < self.y_size:
                        self.grid[x, y].neighbors.append((self.grid[nx, ny], 1))
                    elif self.wrap:
                        if nx < 0:
                            nx = self.x_size - 1
                        elif nx >= self.x_size:
                            nx = 0
                        if ny < 0:
                            ny = self.y_size - 1
                        elif ny >= self.y_size:
                            ny = 0
                        self.grid[x, y].neighbors.append((self.grid[nx, ny], 1))

    def clear_timestamps(self):
        """
        Set the timestamps of all nodes in the grid to 0.
        """
        for x in range(self.x_size):
            for y in range(self.y_size):
                self.grid[x, y].timestamp = 0

    def set_global_state(self, state_grid):
        """
        Set the states of nodes using a 2D array or numpy array of state
        strings. Also resets timestamps.
        """
        if isinstance(state_grid, list):
            state_grid = np.array(state_grid)
        if state_grid.shape != self.grid.shape:
            warn(Warning("State grid set to state with different size than " +
                         "previously set. Changing size."))
            self.grid = np.empty(state_grid.shape, np.dtype(object))
            self.x_size = state_grid.shape[0]
            self.y_size = state_grid.shape[1]
            self.populate_grid()
        for x in range(self.x_size):
            for y in range(self.y_size):
                self.grid[x, y].state = state_grid[x, y]
        self.clear_timestamps()

    def get_global_state(self):
        """
        Get the global state of nodes as a 2D numpy array of strings.
        """
        state_grid = np.empty([self.x_size, self.y_size], np.dtype(object))
        for x in range(self.x_size):
            for y in range(self.y_size):
                state_grid[x, y] = self.grid[x, y].state
        return state_grid

    def getnode(self, x, y):
        return self.grid[x, y]

    def get_node_by_id(self, node_id):
        for node in self:
            if node.node_id == node_id:
                return node

    def species_count(self):
        """
        :return: count of each species in the grid
        """
        return Counter([str(n.state) for n in self])

    @property
    def num_nodes(self):
        return sum([1 for _ in self])

    def __iter__(self):
        return SquareGridIterator(self)

    def __str__(self):
        ret_str = str(self.x_size) + " x " + str(self.y_size) + \
                  " simulation grid:"
        for y in range(self.y_size):
            ret_str += "\n["
            for x in range(self.x_size):
                if x > 0:
                    ret_str += ", "
                ret_str += self.grid[x, y].state
            ret_str += ']'
        return ret_str


class SquareGridIterator:
    def __init__(self, simgrid):
        self.simgrid = simgrid
        self.x = 0
        self.y = 0
        self.done = False

    def __iter__(self):
        return self

    def __next__(self):
        if self.done:
            raise StopIteration
        next_node = self.simgrid.getnode(self.x, self.y)
        self.x += 1
        if self.x >= self.simgrid.x_size:
            self.x = 0
            self.y += 1
            if self.y >= self.simgrid.y_size:
                self.done = True
        return next_node

    def next(self):
        return self.__next__()


class SquareGridWithCornerLeak(SquareGrid):
    """
    A square grid, but where reactions can also happen along corners at some
    (usually small) rate relative to the usual reaction rate.
    """

    def __init__(self, x_size, y_size, corner_rate, wrap=False):
        """
        x_size: The number of nodes from left to right.
        y_size: The number of nodes from top to bottom.
        corner_rate: The speed of reactions along corner connections, as a
                        fraction of the rate of reactions along full edge
                        connections.
        wrap: Iff true, connects the left and right edges together and the top
                and bottom edges together. Default False.
        """
        self.corner_rate = corner_rate
        super(SquareGridWithCornerLeak, self).__init__(x_size, y_size,
                                                       wrap=wrap)

    def populate_grid(self):
        """
        Set/reset all nodes to new nodes with no state and populate neighbor
        nodes. Should only be used by initialization and set routines.
        """
        # Do all the usual connecting...
        super(SquareGridWithCornerLeak, self).populate_grid()

        # ...then add corner connections.
        for x in range(self.x_size):
            for y in range(self.y_size):
                for dx, dy in [(-1, -1), (1, -1), (1, 1), (-1, 1)]:
                    nx = x + dx
                    ny = y + dy
                    if nx >= 0 and nx < self.x_size and \
                            ny >= 0 and ny < self.y_size:
                        neighbor_tuple = (self.grid[nx, ny], self.corner_rate)
                        self.grid[x, y].neighbors.append(neighbor_tuple)
                    elif self.wrap:
                        if nx < 0:
                            nx = self.x_size - 1
                        elif nx >= self.x_size:
                            nx = 0
                        if ny < 0:
                            ny = self.y_size - 1
                        elif ny >= self.y_size:
                            ny = 0
                        neighbor_tuple = (self.grid[nx, ny], self.corner_rate)
                        self.grid[x, y].neighbors.append(neighbor_tuple)


class HexGrid(SquareGrid):
    """
    Representation of a CRN on a 2D hex mesh grid, aligned with hex cell sides
    vertical, and with odd-numbered rows offset to the left and even-numbered
    rows offset to the right.

    Mostly identical to a SquareGrid, except with a different connectivity map.
    Because of this, the overall grid must be "rectangular" -- that is, every
    row must have the same number of cells.
    """

    def __init__(self, x_size, y_size, wrap=False):
        if wrap and y_size % 2 == 1:
            raise ValueError("Can't make a wrapping hex grid with an odd " + \
                             "number of rows. It just doesn't work out.")
        super(HexGrid, self).__init__(x_size, y_size, wrap=wrap)

    def populate_grid(self):
        """
        Set/reset all nodes to new nodes with no state and populate neighbor
        nodes. Should only be used by initialization and set routines.
        """
        for x in range(self.x_size):
            for y in range(self.y_size):
                self.grid[x, y] = (Node(), 1)
                self.grid[x, y].position = (x, y)
        # Populate node neighbor lists
        for x in range(self.x_size):
            for y in range(self.y_size):
                relative_idxs = [(0, -1), (0, 1), (-1, 0), (1, 0)]
                if y % 2 == 1:
                    relative_idxs.append((1, -1))
                    relative_idxs.append((1, 1))
                else:
                    relative_idxs.append((-1, -1))
                    relative_idxs.append((-1, 1))
                for dx, dy in relative_idxs:
                    nx = x + dx
                    ny = y + dy
                    if nx >= 0 and nx < self.x_size and \
                            ny >= 0 and ny < self.y_size:
                        self.grid[x, y].neighbors.append(
                            (self.grid[nx, ny], 1))
                    elif self.wrap:
                        if nx < 0:
                            nx = self.x_size - 1
                        if nx >= self.x_size:
                            nx = 0
                        if ny < 0:
                            ny = self.y_size - 1
                        if ny >= self.y_size:
                            ny = 0
                        self.grid[x, y].neighbors.append(
                            (self.grid[nx, ny], 1))
