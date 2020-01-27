class Substance:
    """
    Defines a substance in a reaction
    """
    def __init__(self, symbol: str, binding_energy: float):
        self.symbol = symbol
        if binding_energy:
            self.binding_energy = binding_energy
        else:
            # TODO: set default binding energy based on symbol
            self.binding_energy = 600
