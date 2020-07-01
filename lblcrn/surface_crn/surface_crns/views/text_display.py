from lblcrn.surface_crn.surface_crns.base import *
import pygame

class TextDisplay(object):
    '''
    Displays text. Encapsulates various formatting decisions for pygame text
    display.

    This object is assigned a width on creation. Text is centered within that
    space.
    '''
    # TODO: make sure this works across all systems
    FONT_NAME         = 'ariaboldttc'  # This is good for MacOS
    BLACK             = (0,0,0)
    WHITE             = (255,255,255)
    HORIZONTAL_BUFFER = 5
    VERTICAL_BUFFER   = 5


    def __init__(self, width, font_size=24, text=""):
        debug = False

        self.display_width = width

        self.font_name = TextDisplay.FONT_NAME
        self.font = pygame.font.SysFont(self.font_name, font_size)

        # Calculate this display's size
        self.display_height = self.font.get_linesize() + 2 * TextDisplay.VERTICAL_BUFFER

        self.text = text


    def update_text(self, text = None):
        if text:
            self.text = text
            return

        font_size = 24
        temp_font = self.font
        while temp_font.size(self.text)[0] > self.display_width:
            font_size *= 0.75
            temp_font = pygame.font.SysFont(self.font_name, font_size)
        self.text_surface = temp_font.render(self.text, True, TextDisplay.BLACK)
        self.text_box = self.text_surface.get_rect()
        self.text_box.centerx = int(self.display_width / 2)
        self.text_box.centery = int(self.display_height / 2)

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

        if debug:
            print("Displaying text at position (" + str(self.x_pos) +
                  ", " + str(self.y_pos) + "), width = " +
                  str(self.display_width) + ", height " +
                  str(self.display_height))

        self.parent_surface = parent_surface
        self.display_surface = self.parent_surface.subsurface(
                          x_pos, y_pos, self.display_width, self.display_height)
        self.display_surface.fill(TextDisplay.WHITE)
        self.display_surface.blit(self.text_surface, self.text_box)

    def get_text(self):
        return self._text

    def set_text(self, value):
        debug = False
        if debug:
            print("Updating display text to " + str(value))
        self._text = value
        self.update_text()

    text = property(get_text, set_text)
# end class TextDisplay
