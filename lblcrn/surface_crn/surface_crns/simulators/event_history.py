import numpy as np

class EventHistory(object):
    '''
    Holds a time-course of surface CRN reaction events. Iterating through an
    EventHistory produces identical results to those of running a simulator
    object, except that the list of events is frozen and deterministic.

    An EventHistory object keeps a pointer to a "current" reaction; this pointer
    can be updated in a forward or a backward direction.
    '''
    def __init__(self):
        self.history = []
        self.idx     = -1

    def at_beginning(self):
        '''
        Returns True iff the index is at the beginning of the event history.
        '''
        return self.idx == -1

    def at_end(self):
        '''
        Returns True iff the index is at the end of the event history.
        '''
        return self.idx + 1 == len(self.history)

    def add_event(self, event):
        '''
        Appends a new event to the end of the eventHistory.
        '''
        self.history.append(event)

    def increment_event(self, distance):
        '''
        Moves the current pointer distance away from the current pointer, then
        returns the event at that pointer.

        Positive values of distance move the pointer forward; negative values
        move the pointer backward.

        Raises an exception if the jump would take the pointer outside the
        eventHistory's history.
        '''
        new_idx = self.idx + distance
        if new_idx < -1 or new_idx >= len(self.history):
            raise IndexError(("eventHistory with %d events at index %d " + \
                              "attempted to jump by %d positions.") % \
                             (len(self.history), self.idx, distance))
        self.idx = new_idx
        return self.history[self.idx]

    def next_event(self):
        '''
        Returns the event after the current pointer.

        Returns None if at the end of the event history.
        '''
        if self.at_end():
            return None
        return self.history[self.idx + 1]

    def previous_event(self):
        '''
        Returns the event before the current pointer and updates the pointer.

        Returns None if at the end of the event history.
        '''
        if self.at_beginning():
            return None
        return self.history[self.idx]

    def clip(self):
        '''
        Deletes all events after the current event, making the current event the
        last event.

        Intended to be used to "unfreeze" a simulation, allowing it to run with
        new random reactions.
        '''
        if not self.at_end():
            del self.history[self.idx+1:]


