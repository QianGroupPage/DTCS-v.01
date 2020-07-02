from lblcrn.surface_crn.surface_crns.base import *
from lblcrn.surface_crn.surface_crns.constants import *
from lblcrn.surface_crn.surface_crns.ordered_dict import OrderedDict
import pygame

class LegendDisplay:
    '''
    Displays a legend of all of the states and colors used in the simulation.
    '''
    pygame.font.init()
    LEGEND_FONT       = pygame.font.SysFont('menlottc', 16)
    HORIZONTAL_BUFFER = 10 #5
    VERTICAL_BUFFER   = 5
    MIN_BOX_HEIGHT    = 10
    LINE_HEIGHT       = max(MIN_BOX_HEIGHT,
                            LEGEND_FONT.get_linesize() + 2 * VERTICAL_BUFFER)

    def __init__(self, colormap, min_x = 0, min_y = 0):
        debug = False

        self.min_x    = min_x
        self.min_y    = min_y

        # Make a deep copy of colormap so it doesn't get modified
        colormap = OrderedDict(colormap)
        # Apply state classes to reduce colormappings
        if COLOR_CLASSES in colormap:
            classmap = colormap[COLOR_CLASSES]
            for state_class in classmap:
                if debug:
                    print("state_class: " + str(state_class))
                    print("\t classmap[state_class]: " +
                          str(classmap[state_class]))
                    print("\t classmap[state_class][0]: " +
                          str(classmap[state_class][0]))
                class_color = colormap[classmap[state_class][0]]
                colormap[state_class] = class_color
                for state in classmap[state_class]:
                    del colormap[state]

        # Reverse the colormap to build a mapping from colors to strings of
        # states.
        self.statemap = OrderedDict()
        for state, color in colormap.iteritems():
            if state is COLOR_CLASSES:
                continue
            if isinstance(color, list):
                color = map(lambda num: "{0:0>2x}".format(int(num)), color)
                color_str = "0x"
                for substr in color:
                    color_str += substr
                color = color_str
            if color in self.statemap.keys():
                self.statemap[color] += ", " + state
            else:
                self.statemap[color] = state

        # Generate text surfaces for each color
        self.text_surfaces = OrderedDict()
        BLACK = (0,0,0)
        for color in self.statemap.keys():
            # If color is a list of values, covert to a hex string so that a)
            # it is hashable and b) pygame can read it.
            state_text = self.statemap[color]
            state_text_surface = LegendDisplay.LEGEND_FONT.render(state_text,
                                                                  True,
                                                                  BLACK)
            self.text_surfaces[color] = state_text_surface

        # Calculate this display's size
        self.box_height = LegendDisplay.LINE_HEIGHT - \
                          2*LegendDisplay.VERTICAL_BUFFER
        self.box_width  = self.box_height
        self.display_height = self.LINE_HEIGHT * \
                              len(self.text_surfaces.keys())
        self.display_width  = 2*LegendDisplay.HORIZONTAL_BUFFER + \
                              self.box_width + \
                              max(map(
                                lambda color:
                                    self.text_surfaces[color].get_rect().width,
                                self.text_surfaces.keys()))
        if debug:
            print()

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
            print("Displaying color legend at position (" + str(x_pos) + "," +
                  str(y_pos) + ") with width " + str(self.display_width) +
                  " and height " + str(self.display_height) + ".")
        self.parent_surface = parent_surface
        self.display_surface = parent_surface.subsurface(
                        (x_pos, y_pos, self.display_width, self.display_height))

        # Initial render
        n = 0
        for color in self.text_surfaces.keys():
            # Render the color box for this color
            if debug:
                print("Rendering color box with color " + color)
            box_color = pygame.Color(color)
            box_x = LegendDisplay.HORIZONTAL_BUFFER
            box_y = n*LegendDisplay.LINE_HEIGHT + LegendDisplay.VERTICAL_BUFFER
            color_rect = pygame.Rect(box_x, box_y,
                                     self.box_width, self.box_height)
            pygame.draw.rect(self.display_surface, pygame.Color(0,0,0),
                             color_rect, 1)
            if debug:
                print("...at position " + str(color_rect.topleft) +
                      " with size " + str(color_rect.width) + "x" +
                      str(color_rect.height))
            color_rect.left += 1
            color_rect.top += 1
            color_rect.width -= 2
            color_rect.height -= 2
            pygame.draw.rect(self.display_surface, box_color, color_rect)

            # Render the text label for this color
            text_surface = self.text_surfaces[color]
            text_rect = text_surface.get_rect()
            text_x    = LegendDisplay.HORIZONTAL_BUFFER*2 + self.box_width
            text_y    = box_y
            text_rect.topleft = (text_x, text_y)
            self.display_surface.blit(text_surface, text_rect)

            n += 1
# end class LegendDisplay