class TransitionRule(object):
    '''
    Represents one CRN transition rule of the form
    [inputs] -> [outputs]
    with some rate (1/(average clock ticks between consecutive reactions)).
    Chemical species are represented as strings.
    Transition rules are static.
    '''
    def __init__(self, inputs = None, outputs = None, rate = 1):
        if inputs == None:
            self.inputs = []
        if outputs == None:
            self.outputs = []
        if not isinstance(inputs, list):
            self.inputs = [inputs]
        else:
            self.inputs = inputs
        if not isinstance(outputs, list):
            self.outputs = [outputs]
        else:
            self.outputs = outputs
        self.rate = rate

    def __str__(self):
        str_rep = ""
        str_rep += self.inputs[0]
        for i in range(1, len(self.inputs)):
            str_rep += " + " + self.inputs[i]
        str_rep += " --> " + self.outputs[0]
        for i in range(1, len(self.outputs)):
            str_rep += " + " + self.outputs[i]
        str_rep += " (rate: " + str(self.rate) + ")"
        return str_rep

    def __repr__(self):
        return "{inputs:" + str(self.inputs) + \
                ";outputs:" + str(self.outputs) + \
                ";rate:" + str(self.rate) + "}"

    def __eq__(self, other_rule):
        '''
        Define two transition rules to be equal if and only if they have the 
        same set of input species and the same set of output species. two
        transition rules do NOT have to have the same rate constants to be 
        equal.
        '''
        return self.inputs == other_rule.inputs and \
                self.outputs == other_rule.outputs

    def __neq__(self, other_rule):
        return not self == other_rule
#end class TransitionRule