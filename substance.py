class Substance:
    """
    Defines a substance in a reaction
    """
    def __init__(self, symbol: str, binding_energy):
        self.symbol = symbol
        self.binding_energy = binding_energy
