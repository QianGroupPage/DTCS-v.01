import math

import pygame


class SpectrumDisplay(object):
    debug = False
    """
    Displays a Spectrum object.
    """

    def __init__(self, grid, colormap, text, min_x=0, min_y=0, width=-1, height=-1):
        """

        :param grid: the SquareGrid object displayed
        :param colormap: Dictionary defining what colors are assigned to each state
        :param text:
        :param min_x, min_y: Minimum x and y sizes, in pixels, of the surface
                            created. If the natural size of the grid isn't big
                            enough, the grid will be centered and whitespace
                            will be added to fill the excess.
        :param width:
        :param height:
        """
        # Constants
        self.grid_buffer = 5

        # Argument inits
        self.grid = grid
        self.colormap = colormap.copy()
        for k, v in colormap.items():
            # TODO (Ye): versify that this is necessary.
            self.colormap[k] = tuple([c / 255 for c in v])
        self.min_x = min_x
        self.min_y = min_y
        self.max_x = 1800  # 1800
        self.max_y = 1200  # 1200

        self.text_display_text = text
        self.fig_dpi = 100
        self.recalculate_display_sizes()

    def recalculate_display_sizes(self):

        # TODO: Try use the following for the aspect ratio of the pics
        # Calculate some internal variables
        self.total_grid_width = min(int(2 * self.grid_buffer + \
                                        (self.grid.x_size + 0.5) * self.node_width), self.max_x)
        # Height = top buffer + bottom buffer + (height of one whole hex) +
        #           (row-height for each except the first row)
        self.total_grid_height = min(int(2 * self.grid_buffer + \
                                         self.node_width / math.cos(math.pi / 6) +
                                         (self.grid.y_size - 1) * self.node_height), self.max_y)
        # self.total_grid_width = 0
        # self.total_grid_height = 0
        # Total height and width must be at least big enough to fit other
        # elements of the UI.
        self.display_width = max(self.total_grid_width, self.min_x)
        self.display_height = max(self.total_grid_height, self.min_y)

    def render(self, parent_surface, x_pos=0, y_pos=0, width=-1, height=-1):
        debug = False
        '''
        Set up the display and make the first render. This must be called before
        any other updates.
            parent_surface: The surface onto which this grid will be displayed.
            x_pos, y_pos: X and Y coordinates of the upper-left corner of this
                            grid relative to parent_surface.
        '''
        self.x_pos = x_pos - 200
        self.y_pos = y_pos

        # print(f"Desired width {width}; desired height {height}")

        # print("Game grid starting position", x_pos, y_pos, f"width={self.display_width}, height={self.display_height}")
        # Create display surface
        if debug:
            print("Displaying grid display at position (" + str(x_pos) + "," +
                  str(y_pos) + ") with width " + str(self.display_width) +
                  " and height " + str(self.display_height) + ".")
        self.parent_surface = parent_surface

        display_width = self.display_width
        max_height_to_width = display_width / width * height

        display_height = self.display_height
        max_width_to_height = display_height / height * width

        # print(f"Grid picture width={max_width_to_height}, height={display_height}")

        actual_width = display_width
        actual_height = max_height_to_width
        self.total_grid_height = actual_height
        self.total_grid_width = actual_width

        # self.display_surface = parent_surface.subsurface((x_pos, y_pos, actual_width, actual_height))
        self.display_surface = parent_surface.subsurface((x_pos, y_pos, self.display_width, self.display_height))

        # Initial render
        self.update_node(None, full_render=True)

    # TODO: fix the bug of no rendering
    def update_node(self, node, full_render=False):
        '''
        Redraw a specified node, in accordance with standard display node API.

        In practice, node is not needed. We use matplotlib to draw everything. When the simulatpr is updating
        a node, don't do anything.

        TODO: research on and improve speed.
        '''
        if not full_render:
            return
        pic_str, size = self._get_voronoi_pic_string()
        voronoi_image = pygame.image.fromstring(pic_str, size, "RGB")
        self.display_surface.blit(voronoi_image, (self.x_pos, self.y_pos))

    def re_render(self):
        self.update_node(None, full_render=True)
