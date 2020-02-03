from typing import List, Tuple

class Substance:
    """
    Defines a substance in a reaction.

    Schedule is a list of 2-length tuples with the first value denoting the time of change,
    and the second denoting the actual amount.
    """
    def __init__(self, symbol: str, final_val: float, schedule: List[Tuple], binding_energy: float):
        self.symbol = symbol
        if binding_energy:
            self.binding_energy = binding_energy
        else:
            # TODO: set default binding energy based on symbol
            self.binding_energy = 600
        self.final_val = final_val
        self.schedule = schedule

