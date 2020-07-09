from lblcrn.surface_crn.surface_crns.base import *
from lblcrn.surface_crn.surface_crns.views.text_display import TextDisplay
import pygame

class TimeDisplay(TextDisplay):
    '''
    Displays a text description of the time.

    This object is assigned a width on creation. Text is centered within that
    space.
    '''
    def __init__(self, width):
        super(TimeDisplay, self).__init__(width, font_size=18)
        self.time = 0

    def update_time_text(self):
        self.text = "T = {0:.2f}".format(self.time)

    def get_time(self):
        return self._time

    def set_time(self, value):
        debug = False
        if debug:
            print("Updating time text time to " + str(value))
        self._time = value
        self.update_time_text()

    time = property(get_time, set_time)
# end class TimeDisplay