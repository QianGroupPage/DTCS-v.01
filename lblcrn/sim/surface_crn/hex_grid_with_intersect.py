import math
import warnings

import numpy as np
import pygame

from lblcrn.sim.surface_crn.surface_crns.base.node import Node
from lblcrn.sim.surface_crn.surface_crns.models.grids import HexGrid
from lblcrn.sim.surface_crn.surface_crns.views.grid_display import HexGridDisplay


class HexGridPlusIntersections(HexGrid):
    """
    Represents a rectangular hex grid (as in HexGrid), plus a site at each
    3-fold intersection site between major sites. Only intersection sites with
    all three node neighbors present are modeled -- intersections along the
    edges are omitted (unless the grid is set to wrap).

    Code for this class borrows heavily from the
    surface_crns.models.grids.SquareGrid and surface_crns.models.grids.HexGrid
    classes. If you want to make a lattice that is very similar to one of those
    two, you may want to consider subclassing either SquareGrid or HexGrid
    instead.

    The SquareGrid and HexGrid classes store their nodes in a rectangular numpy
    array. In addition, this class stores a second (also rectangular numpy array
    of the intersection sites.
    """
    def __init__(self, x_size, y_size, wrap = False):
        if not isinstance(x_size, int) or not isinstance(y_size, int):
            # TODO: use TypeException here
            raise Exception("SimGrid dimensions must be integers")
        if wrap and y_size%2 == 1:
            raise Exception("Can't make a wrapping hex grid with an odd " + \
                            "number of rows. It just doesn't work out.")

        self.x_size     = x_size
        self.y_size     = y_size
        self.int_x_size = (self.x_size - 1)*2  # #columns of intersection sites
        self.int_y_size = self.x_size - 1  # #rows of intersection sites
        self.wrap       = wrap
        self.grid       = np.empty([x_size, y_size], np.dtype(object))
        self.grid       = np.array(self.grid)
        self.int_grid   = np.empty([self.int_x_size, self.int_y_size],
                                    np.dtype(object))
        self.int_grid = np.array(self.int_grid)
        # Instantiate nodes
        self.populate_grid()

    def populate_grid(self):
        """
        Set/reset all nodes to new nodes with no state and populate neighbor
        nodes. Should only be used by initialization and set routines.
        """
        # Initialize both grids with no neighbors (yet).
        for x in range(self.x_size):
            for y in range(self.y_size):
                new_node = Node()
                new_node.position = (x,y)
                new_node.is_intersection = False
                self.grid[x,y] = new_node
        for x in range(self.int_x_size):
            for y in range(self.int_y_size):
                new_node = Node()
                new_node.position = (x, y) # Different x-scale than the hex grid
                new_node.is_intersection = True
                self.int_grid[x,y] = new_node

        # Populate node neighbor lists of hex grid nodes
        for x in range(self.x_size):
            for y in range(self.y_size):
                relative_idxs = [(0, -1), (0,1), (-1, 0), (1, 0)]
                if y%2 == 1:
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
                        self.grid[x,y].neighbors.append((self.grid[nx, ny], 1))
                    elif self.wrap:
                        if nx < 0:
                            nx = self.x_size-1
                        if nx >= self.x_size:
                            nx = 0
                        if ny < 0:
                            ny = self.y_size-1
                        if ny >= self.y_size:
                            ny = 0
                        self.grid[x,y].neighbors.append((self.grid[nx, ny], 1))
        # Populate node neighbor lists of intersection nodes
        for x in range(self.int_x_size):
            for y in range(self.int_y_size):
                # Can classify intersection nodes into two categories: those
                # on the top point of a hex, and those on the bottom point of
                # a hex. The local maps for these two types are rather
                # different, so we'll treat them separately.
                top_point = (x+y)%2 == 0
                if top_point:
                    relative_idxs = [(0, 0),
                                     (1, 0),
                                     (y%2, 1)]
                else:
                    relative_idxs = [((y+1)%2, 0),
                                     (1, 1),
                                     (0, 1)]
                for dx, dy in relative_idxs:
                    nx = int(x/2) + dx
                    ny = y + dy
                    self.grid[nx, ny].neighbors.append((self.int_grid[x, y], 1))
                    self.int_grid[x, y].neighbors.append((self.grid[nx, ny], 1))

    def get_hex_node(self, x, y):
        return self.grid[x,y]

    def get_intersection_node(self, x, y):
        return self.int_grid[x,y]

    def set_global_state(self, state_grid, state_grid_3f=None):
        """
        Set the states of nodes using a 2D array or numpy array of state
        strings. Also resets timestamps.
        """
        state_grid = np.asarray(state_grid)
        state_grid_3f = np.asarray(state_grid_3f)

        if state_grid.shape != self.grid.shape:
            raise ValueError(f'Incorrect top site state grid size '
                             f'{state_grid.shape}, should '
                             f'match .grid.shape ({self.grid.shape})')
        elif state_grid_3f.shape != self.int_grid.shape:
            raise ValueError('Incorrect threefold site state grid size  '
                             f'{state_grid_3f.shape}, should '
                             f'match .grid.shape ({self.int_grid.shape})')

        for x in range(self.x_size):
            for y in range(self.y_size):
                self.grid[x, y].state = state_grid[x, y]
        for x in range(self.int_x_size):
            for y in range(self.int_y_size):
                self.int_grid[x, y].state = state_grid_3f[x, y]

        self.clear_timestamps()

    def __iter__(self):
        return HexPlusIntersectionGridIterator(self)

    def __str__(self):
        ret_str = str(self.x_size) + " x " + str(self.y_size) + \
                    " hex grid:"
        for y in range(self.y_size):
            ret_str += "\n["
            for x in range(self.x_size):
                if x > 0:
                    ret_str += ", "
                ret_str += self.grid[x,y].state
            ret_str += ']\nintersection grid:'
        for y in range(self.int_y_size):
            ret_str += "\n["
            for x in range(self.int_x_size):
                if x > 0:
                    ret_str += ", "
                ret_str += self.int_grid[x,y].state
            ret_str += ']'
        return ret_str


class HexPlusIntersectionGridIterator:
    """
    Iterate through hex grid, then intersection grid.
    """
    def __init__(self, simgrid):
        self.simgrid = simgrid
        self.x = 0
        self.y = 0
        self.hex_done          = False
        self.intersection_done = False

    def __iter__(self):
        return self

    def __next__(self):
        if self.intersection_done:
            raise StopIteration

        if self.hex_done:
            next_node = self.simgrid.get_intersection_node(self.x, self.y)
        else:
            next_node = self.simgrid.getnode(self.x, self.y)
        self.x += 1
        if (not self.hex_done and self.x >= self.simgrid.x_size) or \
           (self.hex_done and self.x >= self.simgrid.int_x_size):
            self.x = 0
            self.y += 1
            if (not self.hex_done and self.y >= self.simgrid.y_size) or \
               (self.hex_done and self.y >= self.simgrid.int_y_size):
                if not self.hex_done:
                    self.hex_done = True
                    self.x = 0
                    self.y = 0
                else:
                    self.intersection_done = True
        return next_node

    def next(self):
        return self.__next__()


class HexGridPlusIntersectionDisplay(HexGridDisplay):
    """
    See documentation in surface_crns.views.grid_display.SurfaceDisplay for
    details on how to make a valid display class. Briefly, a SurfaceDisplay
    needs to override the __init__, render, and update_node functions.

    This particular class subclasses the HexGridDisplay clas (in the
    same file as SurfaceDisplay; a lot of the display code is re-used from that.
    In brief, this class uses the same rendering function as HexGridDisplay,
    then draws another layer on top of that with the intersection nodes. The
    only change between this class and HexGridDisplay is that this class has
    an extended render function that knows how to handle intersection nodes.
    """
    def update_node(self, node):
        """
        This function should blit a single node, or otherwise do whatever is
        required when a node changes state. This will be called whenever the
        state of a node is changed, passing the node that changed.

        params:
            node: The node that just changed state.
        """
        if node.is_intersection:
            center = self.get_center(node)
            center = list(map(int, center))
            radius = int(self.pixels_per_node[0] / 4)
            node_color = self.colormap[node.state]
            pygame.draw.circle(self.display_surface, node_color, center, radius)
            pygame.draw.circle(self.display_surface, (0,0,0), center, radius,
                               1)
            if self.display_text:
                node_text_surface = self.make_node_text(node)
                text_rect = node_text_surface.get_rect()
                text_rect.center = center
                self.display_surface.blit(node_text_surface,
                                          text_rect)
        else:
            super(HexGridPlusIntersectionDisplay, self).update_node(node)
            for n, weight in node.neighbors:
                if n.is_intersection:
                    self.update_node(n)
            center = self.get_center(node)

    def get_center(self, node):
        """
        Returns the coordinates (in pixesls) of the center of this node.
        """
        if not node.is_intersection:
            return super(HexGridPlusIntersectionDisplay, self).get_center(node)

        # As with grid setup, will split this into two cases -- one where the
        # intersection is at the top corner of a hex, and one where it is at the
        # bottom corner of a hex.
        node_x = node.position[0]
        node_y = node.position[1]
        center_to_top_dist = self.node_width * 0.5 / math.cos(math.pi/6.0)
        is_top_corner = (node_x + node_y)%2 == 0
        if is_top_corner:
            center_below = self.get_center(self.grid.get_hex_node(
                                                    int(node_x/2) + node_y%2,
                                                    node_y+1))
            x = center_below[0]
            y = center_below[1] - center_to_top_dist
        else:

            center_above = self.get_center(self.grid.get_hex_node(
                                                int(node_x/2) + (node_x)%2,
                                                node_y))
            x = center_above[0]
            y = center_above[1] + center_to_top_dist

        return (x,y)