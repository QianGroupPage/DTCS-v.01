class Event:
    def __init__(self, **kwargs):
        self.time = kwargs['time'] # Time that the reaction occurred.
        self.rule = kwargs['rule']
        self.participants  = kwargs['participants']
        self.time_issued   = kwargs['time_issued'] # Time when the simulator
                                                   # planned this event.

    def __str__(self):
        ret_str = "Event{"
        ret_str += "time:" + str(self.time) + "; "
        ret_str += "rule:" + str(self.rule) + "; "
        ret_str += "time_issued:" + str(self.time_issued) + "; "
        ret_str += "participants:" + str(self.participants) + "}"
        return ret_str

    def __lt__(self, other):
        return self.time < other.time