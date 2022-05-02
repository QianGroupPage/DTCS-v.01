import math

import matplotlib
import matplotlib.backends.backend_agg as agg
import pygame


class CoordGridDisplay(object):
    debug = False
    """
    Displays a HexGrid object as a colored honeycomb.
    """

    def __init__(self, grid, colormap, min_x=0, min_y=0, width=-1, height=-1, pixels_per_node=5,
                 display_text=False):
        """
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
        self.pixels_per_node = pixels_per_node

        self.display_text = display_text

        self.fig_dpi = 100

        self.recalculate_display_sizes()
        # print(self.display_width, self.display_height)

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

    def pixels_per_node():
        doc = "The width, in pixels, of a single hex in the grid. " + \
              "Setting this also sets the row height (the number of vertical " + \
              "pixels added by adding a row)."

        def fget(self):
            return (self.node_width, self.node_height)

        def fset(self, value):
            if isinstance(value, int):
                self.node_width = value
                self.node_height = value / 2 / math.tan(math.pi / 6)
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

    # TODO: render coord_grid.voronoi_pic
    def render(self, parent_surface, x_pos=0, y_pos=0, width=-1, height=-1):
        debug = False
        """
        Set up the display and make the first render. This must be called before
        any other updates.
            parent_surface: The surface onto which this grid will be displayed.
            x_pos, y_pos: X and Y coordinates of the upper-left corner of this
                            grid relative to parent_surface.
        """
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

    @property
    def _grid_size(self):
        return (self.total_grid_width, self.total_grid_height)

    def _get_voronoi_pic_string(self):
        backend = matplotlib.rcParams['backend']
        matplotlib.use("Agg")
        # TODO: determine if this solves the issue with inconsistent font.
        matplotlib.rcParams["font.family"] = "arial"
        import matplotlib.pyplot as plt
        # Set the fig_size to be the same as the display surface.
        plt.gcf().set_dpi(self.fig_dpi)
        plt.gcf().set_figheight((self.total_grid_height - 2 * self.grid_buffer) / self.fig_dpi)
        plt.gcf().set_figwidth((self.total_grid_width - 2 * self.grid_buffer) / self.fig_dpi)

        # Block to prevent white spaces from showing
        plt.gca().set_axis_off()
        plt.subplots_adjust(top=1, bottom=0, right=1, left=0,
                            hspace=0, wspace=0)
        plt.margins(0, 0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        # NullLocator    plt.savefig("filename.pdf", bbox_inches='tight',
        #                 pad_inches=0)
        #     print("intended", self.total_grid_width, self.total_grid_height)
        # TODO: fix the sizing
        fig = self.grid.voronoi_pic(ax=plt.gca(), return_fig=True, color_index=self.colormap, set_fig_size=False)
        # fig.tight_layout()
        # TODO: currently it seems looping to create this.
        # print("fig size", fig.get_size_inches()*fig.dpi)
        # matplotlib.pyplot.ylim((0, y_upper_limit))
        canvas = agg.FigureCanvasAgg(fig)
        canvas.draw()
        renderer = canvas.get_renderer()
        matplotlib.use(backend)
        # TODO: close plt in case too many plots got opened.
        return renderer.tostring_rgb(), canvas.get_width_height()

    # TODO: fix the bug of no rendering
    def update_node(self, node, full_render=False):
        """
        Redraw a specified node, in accordance with standard display node API.

        In practice, node is not needed. We use matplotlib to draw everything. When the simulatpr is updating
        a node, don't do anything.

        TODO: research on and improve speed.
        """
        if not full_render:
            return
        pic_str, size = self._get_voronoi_pic_string()
        voronoi_image = pygame.image.fromstring(pic_str, size, "RGB")
        self.display_surface.blit(voronoi_image, (self.x_pos, self.y_pos))

    def re_render(self):
        self.update_node(None, full_render=True)

# def make_node_hex(self, node):
#         """
#         Returns the list of vertices of the hex at the node's position.
#         """
#         debug = False
#
#         x_pos, y_pos = self.get_center(node)
#         a = self.node_width * 0.5 / math.cos(math.pi/6.0)
#         b = self.node_width * 0.5 * math.tan(math.pi/6.0)
#         vertex_list = [(x_pos, y_pos + a),
#                        (x_pos + 0.5 * self.node_width, y_pos + b),
#                        (x_pos + 0.5 * self.node_width, y_pos - b),
#                        (x_pos, y_pos - a),
#                        (x_pos - 0.5 * self.node_width, y_pos - b),
#                        (x_pos - 0.5 * self.node_width, y_pos + b)]
#         vertex_list = list(map(lambda pair: (int(pair[0]), int(pair[1])),
#                                vertex_list))
#         if self.debug:
#             print("Making new polygon (hex) with the following vertices: " + \
#                   str(vertex_list))
#
#         return vertex_list
#
#     def get_center(self, node):
#         """
#         Returns the coordinates (in pixesls) of the center of this node.
#         """
#         x = node.position[0]
#         y = node.position[1]
#         # Grid might be floating in a space required by other UI elements.
#         # If so, add a buffer to each side.
#         if self.total_grid_width < self.min_x:
#             x_buffer = (self.min_x - self.total_grid_width)/2
#         else:
#             x_buffer = 0
#         if self.total_grid_height < self.min_y:
#             y_buffer = (self.min_y - self.total_grid_height)/2
#         else:
#             y_buffer = 0
#
#         x_pos = (x + 0.5*(y%2) + 0.5) * self.node_width
#         y_pos = self.node_width * math.tan(math.pi/6.0) + \
#                     y * self.node_width / 2.0 / math.tan(math.pi/6.0)
#         if self.debug:
#             print("Calculated center of node (%d, %d) at (%d, %d)." % \
#                  (x, y, x_pos, y_pos))
#         return (x_pos, y_pos)
#
#     def make_node_text(self, node):
#         BLACK = (0,0,0)
#         WHITE = (255,255,255)
#         node_color = self.colormap[node.state]
#         if sum(node_color) < 150:
#             text_color = WHITE
#         else:
#             text_color = BLACK
#         font = pygame.font.SysFont('monospace', 10)
#         text_surface = font.render(node.state, True, text_color)
#         return text_surface
# end class HexGridDisplay
