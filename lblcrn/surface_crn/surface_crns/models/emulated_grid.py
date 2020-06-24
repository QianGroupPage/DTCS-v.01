import numpy as np
from lblcrn.surface_crn.surface_crns.base.node import *
from lblcrn.surface_crn.surface_crns.models.grids import SquareGrid
from warnings import *

class EmulatedSimulationGrid(object):
    '''
    Representation of a CRN on a 2D finite rectangular mesh grid.
    '''
    def __init__(self, x_size, y_size):
        if not isinstance(x_size, int) or not isinstance(y_size, int):
            raise TypeException("SimGrid dimensions must be integers")
        self.x_size = x_size
        self.y_size = y_size
        self.grid   = np.empty([x_size, y_size], np.dtype(object))
        self.grid   = np.array(self.grid)
        # Instantiate nodes
        self.populate_grid()

    def populate_grid(self):
        '''
        Set/reset all nodes to new nodes with no state and populate neighbor
        nodes. Should only be used by initialization and set routines.
        '''
        for x in range(self.x_size):
            for y in range(self.y_size):
                self.grid[x,y] = Node()
                self.grid[x,y].position = (x,y)
        # Populate node neighbor lists
        for x in range(self.x_size):
            for y in range(self.y_size):
                for dx, dy in [(-1,0), (1,0), (0,-1), (0,1)]:
                    nx = x + dx
                    ny = y + dy
                    if nx >= 0 and nx < self.x_size and \
                        ny >= 0 and ny < self.y_size:
                        self.grid[x,y].neighbors.append(self.grid[nx, ny])

    def clear_timestamps(self):
        '''
        Set the timestamps of all nodes in the grid to 0.
        '''
        for x in range(self.x_size):
            for y in range(self.y_size):
                self.grid[x,y].timestamp = 0

    def set_global_state(self, state_grid):
        '''
        Set the states of nodes using a 2D array or numpy array of state
        strings. Also resets timestamps.
        '''
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
                self.grid[x,y].state = state_grid[x,y]
        self.clear_timestamps()

    def get_global_state(self):
        '''
        Get the global state of nodes as a 2D numpy array of strings.
        '''
        state_grid = np.empty([self.x_size, self.y_size], np.dtype(object))
        for x in range(self.x_size):
            for y in range(self.y_size):
                state_grid[x,y] = self.grid[x,y].state
        return state_grid

    def getnode(self, x, y):
        return self.grid[x,y]

    def __iter__(self):
        return SimulationGridIterator(self)

    def __str__(self):
        ret_str = str(self.x_size) + " x " + str(self.y_size) + \
                    " simulation grid:"
        for y in range(self.y_size):
            ret_str += "\n["
            for x in range(self.x_size):
                if x > 0:
                    ret_str += ", "
                ret_str += self.grid[x,y].state
            ret_str += ']'
        return ret_str

class SimulationGridIterator:
    def __init__(self, simgrid):
        self.simgrid = simgrid
        self.x = 0
        self.y = 0
        self.done = False

    def __iter__(self):
        return self

    def next(self):
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
