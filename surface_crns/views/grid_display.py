from surface_crns.base import *
import pygame
import math
'''
Classes for creating and updating pygame surfaces based on simulation data.
'''

class SurfaceDisplay(object):
    '''
    Interface for displays of surface CRNs. If you want to make your own custom
    surface display for the simulator, inherit from this class and override all
    of its methods.

    Each Surface object also needs to have the following accessible variables:
        * display_width: The total width of the display of the surface, in
                            pixels. For a square grid, this is
                            max(min_x, <# columns> * <cell width in pixels>)
        * display_width: The total height of the display of the surface, in
                            pixels. For a square grid, this is
                            max(min_y, <# columns> * <cell height in pixels>)
        * grid: The surface being simulated. Doesn't actually need to be a
                *grid*, just needs to be iterable (and probably should have
                *some* kind of structure).

    Note that the window in which you will be embedding your display requires at
    least some space to display other elements. These will be given to your
    object as the min_x and min_y parameters in the constructor. Make sure your
    display is at least this big, or you'll get some weird display issues.
    '''
    def __init__(self, grid, colormap, min_x = 0, min_y = 0, pixels_per_node=5,
                 display_text = False):
        '''
        Params:
            grid: Whatever object represents the nodes on the surface. For a
                    square grid, this would be a SquareGrid object.
            colormap: A dictionary defining what color each state will be
                        displayed as. Maps from states (strings) to colors (RGB
                        3-tuples, or anything else directly recognizable by
                        Pygame).
            min_x, min_y: Minimum x and y sizes, in pixels, of the
                            surface-displaying portion of the window. The window
                            requires at least some space to blit control buttons
                            and a legend; make sure the display_width and
                            display_height you calcuate are at least this big,
                            or you will have display problems.
            pixels_per_node: Scaling factor. Usually the size of a single node,
                                in pixels.
            display_text: If true, your display should blit the state of each
                            node xonto that node.
        '''
        raise NotImplementedError("You need to override the constructor for " +\
                                  "Surface.")

    def render(self, parent_surface, x_pos = 0, y_pos = 0):
        '''
        This function should blit the entire surface onto its parent. This will
        be called once at the beginning of the simulation.

        Params:
            parent_surface: The pygame.Surface object representing the entire
                            window.
            x_pos, y_pos: X and Y coordinates of the upper-left corner of this
                            grid relative to parent_surface. Use these to make
                            sure you aren't blitting onto stuff outside the
                            display region.
        '''
        raise NotImplementedError("You need to override the 'render' " + \
                                  "method of Surface.")

    def update_node(self, node):
        '''
        This function should blit a single node, or otherwise do whatever is
        required when a node changes state. This will be called whenever the
        state of a node is changed, passing the node that changed.

        params:
            node: The node that just changed state.
        '''
        raise NotImplementedError("You need to override the 'update_node' " + \
                                  "method of Surface.")


class SquareGridDisplay(object):
    debug = False
    '''
    Displays a SquareGrid object as a colored grid.
    '''

    def __init__(self, grid, colormap, min_x = 0, min_y = 0, pixels_per_node=5,
                 display_text = False):
        '''
         Parameters:
            grid: The SquareGrid object displayed
            colormap: Dictionary defining what colors are assigned to each state
            min_x, min_y: Minimum x and y sizes, in pixels, of the surface
                            created. If the natural size of the grid isn't big
                            enough, the grid will be centered and whitespace
                            will be added to fill the excess.
            pixels_per_node: Width and height of each node, in pixels, either as
                                a pair of integers (for arbitrary form factors)
                                or as a single integer (for square nodes)
            display_text: If True, will display each node with a text overlay of
                            that node's state. Otherwise, will only display the
                            color of the node.
        '''

        # Constants
        self.grid_buffer = 5

        # Argument inits
        self.grid = grid
        self.colormap = colormap
        self.min_x = min_x
        self.min_y = min_y
        self.pixels_per_node = pixels_per_node

        self.display_text = display_text

        self.recalculate_display_sizes()


    def recalculate_display_sizes(self):
        # Calculate some internal variables
        self.display_width  = max((2 * self.grid_buffer +
                          self.grid.x_size * self.node_width), self.min_x)
        self.display_height = max((2 * self.grid_buffer +
                          self.grid.y_size * self.node_height), self.min_y)

    def pixels_per_node():
        doc = "The number of pixels used to display a single node in the grid."
        def fget(self):
            return (self.node_width, self.node_height)
        def fset(self, value):
            if isinstance(value, int):
                self.node_width = value
                self.node_height = value
            elif (isinstance(value, list) or isinstance(value, tuple)) and \
                 len(value) == 2:
                self.node_width = value[0]
                self.node_height = value[1]
            else:
                raise Exception("Invalid argument for pixels_per_node: " +
                                str(value))
            self.recalculate_display_sizes()
        def fdel(self):
            del self.node_width
            del self.node_height
        return locals()
    pixels_per_node = property(**pixels_per_node())

    def render(self, parent_surface, x_pos = 0, y_pos = 0):
        debug = False
        '''
        Set up the display and make the first render. This must be called before
        any other updates.
            parent_surface: The surface onto which this grid will be displayed.
            x_p os, y_pos: X and Y coordinates of the upper-left corner of this
                            grid relative to parent_surface.
        '''
        self.x_pos = x_pos
        self.y_pos = y_pos
        # Create display surface
        if debug:
            print("Displaying grid display at position (" + str(x_pos) + "," +
                  str(y_pos) + ") with width " + str(self.display_width) +
                  " and height " + str(self.display_height) + ".")
        self.parent_surface = parent_surface
        self.display_surface = parent_surface.subsurface(
                        (x_pos, y_pos, self.display_width, self.display_height))
        # Initial render
        for node in self.grid:
            self.update_node(node)

    def update_node(self, node):
        '''
        Redraw a specified node.
        '''
        new_rect   = self.make_node_rectangle(node)
        node_color = self.colormap[node.state]
        pygame.draw.rect(self.display_surface, node_color, new_rect)
        if self.display_text:
            node_text_surface = self.make_node_text(node)
            text_rect = node_text_surface.get_rect()
            text_rect.center = new_rect.center
            self.display_surface.blit(node_text_surface,
                                      text_rect)

    def make_node_rectangle(self, node):
        debug = False

        x = node.position[0]
        y = node.position[1]
        if self.grid.x_size * self.node_width < self.min_x:
            x_buffer = (self.min_x - self.grid.x_size*self.node_width)/2
            x_pos    = x_buffer + x * self.node_width
        else:
            x_pos    = self.grid_buffer + x * self.node_width
        if self.grid.y_size * self.node_height < self.min_y:
            y_buffer = (self.min_y - self.grid.y_size*self.node_height)/2
            y_pos    = y_buffer + y * self.node_height
        else:
            y_pos    = self.grid_buffer + y * self.node_height
        if debug:
            print("Creating new rectangle with position (" + str(x_pos) + "," +
                  str(y_pos) + "), height " + str(self.node_height) +
                  ", and width " + str(self.node_width) + ".")
        return pygame.Rect(x_pos, y_pos, self.node_width, self.node_height)

    def make_node_text(self, node):
        BLACK = (0,0,0)
        WHITE = (255,255,255)
        node_color = self.colormap[node.state]
        if sum(node_color) < 150:
            text_color = WHITE
        else:
            text_color = BLACK
        font = pygame.font.SysFont('monospace', 10)
        text_surface = font.render(node.state, True, text_color)
        return text_surface
#end class SquareGridDisplay


class HexGridDisplay(object):
    debug = False
    '''
    Displays a HexGrid object as a colored honeycomb.
    '''

    def __init__(self, grid, colormap, min_x = 0, min_y = 0, pixels_per_node=5,
                 display_text = False):
        '''
         Parameters:
            grid: The SquareGrid object displayed
            colormap: Dictionary defining what colors are assigned to each state
            min_x, min_y: Minimum x and y sizes, in pixels, of the surface
                            created. If the natural size of the grid isn't big
                            enough, the grid will be centered and whitespace
                            will be added to fill the excess.
            pixels_per_node: Width and height of each node, in pixels, either as
                                a pair of integers (for arbitrary form factors)
                                or as a single integer (for square nodes)
            display_text: If True, will display each node with a text overlay of
                            that node's state. Otherwise, will only display the
                            color of the node.
        '''

        # Constants
        self.grid_buffer = 5

        # Argument inits
        self.grid = grid
        self.colormap = colormap
        self.min_x = min_x
        self.min_y = min_y
        self.pixels_per_node = pixels_per_node

        self.display_text = display_text

        self.recalculate_display_sizes()

    def recalculate_display_sizes(self):
        # Calculate some internal variables
        self.total_grid_width  = int(2 * self.grid_buffer + \
                          (self.grid.x_size+0.5) * self.node_width)
        # Height = top buffer + bottom buffer + (height of one whole hex) +
        #           (row-height for each except the first row)
        self.total_grid_height = int(2 * self.grid_buffer + \
                          self.node_width / math.cos(math.pi/6) +
                          (self.grid.y_size-1) * self.node_height)
        # Total height and width must be at least big enough to fit other
        # elements of the UI.
        self.display_width  = max(self.total_grid_width, self.min_x)
        self.display_height = max(self.total_grid_height, self.min_y)

    def pixels_per_node():
        doc = "The width, in pixels, of a single hex in the grid. " + \
              "Setting this also sets the row height (the number of vertical "+\
              "pixels added by adding a row)."
        def fget(self):
            return (self.node_width, self.node_height)
        def fset(self, value):
            if isinstance(value, int):
                self.node_width = value
                self.node_height = value / 2 / math.tan(math.pi/6)
            elif (isinstance(value, list) or isinstance(value, tuple)) and \
                 len(value) == 2:
                self.node_width = value[0]
                self.node_height = value[1]
            else:
                raise Exception("Invalid argument for pixels_per_node: " +
                                str(value))
            self.recalculate_display_sizes()
        def fdel(self):
            del self.node_width
            del self.node_height
        return locals()
    pixels_per_node = property(**pixels_per_node())

    def render(self, parent_surface, x_pos = 0, y_pos = 0):
        debug = False
        '''
        Set up the display and make the first render. This must be called before
        any other updates.
            parent_surface: The surface onto which this grid will be displayed.
            x_p os, y_pos: X and Y coordinates of the upper-left corner of this
                            grid relative to parent_surface.
        '''
        self.x_pos = x_pos
        self.y_pos = y_pos
        # Create display surface
        if debug:
            print("Displaying grid display at position (" + str(x_pos) + "," +
                  str(y_pos) + ") with width " + str(self.display_width) +
                  " and height " + str(self.display_height) + ".")
        self.parent_surface = parent_surface
        self.display_surface = parent_surface.subsurface(
                        (x_pos, y_pos, self.display_width, self.display_height))
        # Initial render
        for node in self.grid:
            self.update_node(node)

    def update_node(self, node):
        '''
        Redraw a specified node.
        '''
        new_hex    = self.make_node_hex(node)
        node_color = self.colormap[node.state]
        pygame.draw.polygon(self.display_surface, node_color, new_hex)
        pygame.draw.lines(self.display_surface, (0,0,0), True, new_hex)
        if self.display_text:
            node_text_surface = self.make_node_text(node)
            text_rect = node_text_surface.get_rect()
            text_rect.center = self.get_center(node)
            self.display_surface.blit(node_text_surface,
                                      text_rect)

    def make_node_hex(self, node):
        '''
        Returns the list of vertices of the hex at the node's position.
        '''
        debug = False

        x_pos, y_pos = self.get_center(node)
        a = self.node_width * 0.5 / math.cos(math.pi/6.0)
        b = self.node_width * 0.5 * math.tan(math.pi/6.0)
        vertex_list = [(x_pos, y_pos + a),
                       (x_pos + 0.5 * self.node_width, y_pos + b),
                       (x_pos + 0.5 * self.node_width, y_pos - b),
                       (x_pos, y_pos - a),
                       (x_pos - 0.5 * self.node_width, y_pos - b),
                       (x_pos - 0.5 * self.node_width, y_pos + b)]
        vertex_list = list(map(lambda pair: (int(pair[0]), int(pair[1])),
                               vertex_list))
        if self.debug:
            print("Making new polygon (hex) with the following vertices: " + \
                  str(vertex_list))

        return vertex_list

    def get_center(self, node):
        '''
        Returns the coordinates (in pixesls) of the center of this node.
        '''
        x = node.position[0]
        y = node.position[1]
        # Grid might be floating in a space required by other UI elements.
        # If so, add a buffer to each side.
        if self.total_grid_width < self.min_x:
            x_buffer = (self.min_x - self.total_grid_width)/2
        else:
            x_buffer = 0
        if self.total_grid_height < self.min_y:
            y_buffer = (self.min_y - self.total_grid_height)/2
        else:
            y_buffer = 0

        x_pos = (x + 0.5*(y%2) + 0.5) * self.node_width
        y_pos = self.node_width * math.tan(math.pi/6.0) + \
                    y * self.node_width / 2.0 / math.tan(math.pi/6.0)
        if self.debug:
            print("Calculated center of node (%d, %d) at (%d, %d)." % \
                 (x, y, x_pos, y_pos))
        return (x_pos, y_pos)

    def make_node_text(self, node):
        BLACK = (0,0,0)
        WHITE = (255,255,255)
        node_color = self.colormap[node.state]
        if sum(node_color) < 150:
            text_color = WHITE
        else:
            text_color = BLACK
        font = pygame.font.SysFont('monospace', 10)
        text_surface = font.render(node.state, True, text_color)
        return text_surface
#end class HexGridDisplay


class ParallelEmulatedSquareGridDisplay(object):
    debug = False
    '''
    Displays an underlying grid and the grid it emulates (i.e., the underlying
    process of a game of life automaton and the automoton it emulates)
    side-by-side.

    Assumes that each emulated cell is determined by the value of a fixed
    position in the underlying grid; that the first character of the
    representative cell's state gives the emulated cell's state; and that the
    emulated cell's state is constant until the representative cell gains a
    new state beginning with anything other than "B".
    '''

    def __init__(self, grid, colormap, emulation_colormap, horizontal_buffer,
                 vertical_buffer, cell_height, cell_width,
                 representative_cell_x, representative_cell_y, min_x = 0,
                 min_y = 0, pixels_per_node=5, display_text = False):
        '''
         Parameters:
            grid: The SquareGrid object displayed
            colormap: Dictionary defining what colors are assigned to each state
            emulation_colormap: Same as colormap, but for the emulated grid.
            horizontal_buffer, vertical_buffer:
                        The number of underlying grid nodes on the left/right
                        and top/bottom of the underlying grid that aren't
                        included in the emulation (boundary conditions)
            cell_width, cell_height: The number of underlying grid nodes
                                        in an emulated node in each direction.
            representative_cell_x, representative_cell_y:
                        The location, within each emulated cell, of the
                        underlying grid cell containing information about the
                        emulated cell.
            min_x, min_y: Minimum x and y sizes, in pixels, of the surface
                            created. If the natural size of the grid isn't big
                            enough, the grid will be centered and whitespace
                            will be added to fill the excess.
            pixels_per_node: Width and height of each node, in pixels, either as
                                a pair of integers (for arbitrary form factors)
                                or as a single integer (for square nodes)
            display_text: If True, will display each node with a text overlay of
                            that node's state. Otherwise, will only display the
                            color of the node.
        '''
        pygame.init()

        # Constants
        self.grid_buffer = 5


        # Argument inits
        self.grid = grid
        self.horizontal_buffer = horizontal_buffer
        self.vertical_buffer = vertical_buffer
        self.emulated_cell_height = cell_height
        self.emulated_cell_width = cell_width
        self.representative_cell_x = representative_cell_x
        self.representative_cell_y = representative_cell_y
        self.emulation_colormap = emulation_colormap
        self.colormap = colormap
        self.min_x = min_x
        self.min_y = min_y
        self.pixels_per_node = pixels_per_node

        self.display_text = display_text

        self.recalculate_display_sizes()


    def recalculate_display_sizes(self):
        # Calculate some internal variables
        self.display_width  = max((self.grid_buffer +
                                  2 * self.grid.x_size * self.node_width ),
                                 self.min_x)
        self.display_height = max((2 * self.grid_buffer +
                          self.grid.y_size * self.node_height), self.min_y)

    def pixels_per_node():
        doc = "The number of pixels used to display a single node in the grid."
        def fget(self):
            return (self.node_width, self.node_height)
        def fset(self, value):
            if isinstance(value, int):
                self.node_width = value
                self.node_height = value
            elif (isinstance(value, list) or isinstance(value, tuple)) and \
                 len(value) == 2:
                self.node_width = value[0]
                self.node_height = value[1]
            else:
                raise Exception("Invalid argument for pixels_per_node: " +
                                str(value))
            self.recalculate_display_sizes()
        def fdel(self):
            del self.node_width
            del self.node_height
        return locals()
    pixels_per_node = property(**pixels_per_node())

    def render(self, parent_surface, x_pos = 0, y_pos = 0):
        debug = False
        '''
        Set up the display and make the first render. This must be called before
        any other updates.
            parent_surface: The surface onto which this grid will be displayed.
            x_p os, y_pos: X and Y coordinates of the upper-left corner of this
                            grid relative to parent_surface.
        '''
        self.x_pos = x_pos
        self.y_pos = y_pos
        # Create display surface
        if debug:
            print("Displaying grid display at position (" + str(x_pos) + "," +
                  str(y_pos) + ") with width " + str(self.display_width) +
                  " and height " + str(self.display_height) + ".")
        self.parent_surface = parent_surface
        self.display_surface = parent_surface.subsurface(
                        (x_pos, y_pos, self.display_width, self.display_height))
        # Initial render
        for x in range(self.grid.x_size):
            for y in range(self.grid.y_size):
                self.update_node_at_position(x, y)

    def update_node_at_position(self, x, y):
        self.update_node(self.grid.getnode(x,y))

    def update_node(self, node):
        '''
        Redraw a specified node, and its emulated node if that emulated node
        changed state
        '''
        new_rect   = self.make_node_rectangle(node)
        node_color = self.colormap[node.state]
        pygame.draw.rect(self.display_surface, node_color, new_rect)
        if self.display_text:
            node_text_surface = self.make_node_text(node)
            text_rect = node_text_surface.get_rect()
            text_rect.center = new_rect.center
            self.display_surface.blit(node_text_surface,
                                      text_rect)

        # Update the emulated node if necessary
        x = (node.position[0]-self.horizontal_buffer) % self.emulated_cell_width
        y = (node.position[1]-self.vertical_buffer) % self.emulated_cell_height
        if x == self.representative_cell_x and y == self.representative_cell_y:
            represented_x = (node.position[0]-self.horizontal_buffer)          \
                            / self.emulated_cell_width
            represented_y = (node.position[1]-self.vertical_buffer)            \
                            / self.emulated_cell_height
            new_rect = self.make_emulated_node_rectangle(represented_x,
                                                          represented_y,
                                                          node.state[0])
            if node.state[0] in self.emulation_colormap:
                node_color = self.emulation_colormap[node.state[0]]
            else:
                node_color = self.emulation_colormap['B']
            pygame.draw.rect(self.display_surface, node_color, new_rect)

    def make_node_rectangle(self, node):
        debug = False

        x = node.position[0]
        y = node.position[1]
        if self.display_width < self.min_x:
            x_buffer = (self.min_x - self.grid.x_size*self.node_width*2 + \
                                                    self.horizontal_buffer)/2
            x_pos    = x_buffer + x * self.node_width
        else:
            x_pos    = self.grid_buffer + x * self.node_width
        if self.display_height < self.min_y:
            y_buffer = (self.min_y-self.grid.y_size*self.node_height)/2
            y_pos    = y_buffer + y * self.node_height
        else:
            y_pos    = self.grid_buffer + y * self.node_height
        if debug:
            print("Creating new rectangle with position (" + str(x_pos) + "," +
                  str(y_pos) + "), height " + str(self.node_height) +
                  ", and width " + str(self.node_width) + ".")
        return pygame.Rect(x_pos, y_pos, self.node_width, self.node_height)

    def make_emulated_node_rectangle(self, x, y, state):
        if self.display_width < self.min_x:
            x_buffer = (self.min_x - self.grid.x_size*self.node_width*2 + \
                                                    self.horizontal_buffer)/2
        else:
            x_buffer = self.grid_buffer
        x_pos = x_buffer + self.grid.x_size*self.node_width + \
                self.node_width*self.horizontal_buffer + \
                x*self.node_width*self.emulated_cell_width
        if self.display_height < self.min_y:
            y_buffer = (self.min_y-self.grid.y_size*self.node_height)/2
        else:
            y_buffer = self.grid_buffer
        y_pos = y_buffer + (self.node_height * self.vertical_buffer) + \
                (y * self.node_height * self.emulated_cell_height)
        return pygame.Rect(x_pos, y_pos,
                           self.node_width*self.emulated_cell_width,
                           self.node_height*self.emulated_cell_height)

    def make_node_text(self, node):
        BLACK = (0,0,0)
        WHITE = (255,255,255)
        node_color = self.colormap[node.state]
        if sum(node_color) < 150:
            text_color = WHITE
        else:
            text_color = BLACK
        font = pygame.font.SysFont('monospace', 10)
        text_surface = font.render(node.state, True, text_color)
        return text_surface

#end class ParallelEmulatedSquareGridDisplay