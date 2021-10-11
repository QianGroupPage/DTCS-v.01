import monty.json
import sympy as sym


class Marker(monty.json.MSONable):
    """A marker object to mark occurences of a given species with a different name.
    """
    def __init__(self, species: str, name: str, species_symbol: sym.Symbol = None, color: str = ""):
        self.species = species

        self.species_symbol = species_symbol if species_symbol else sym.Symbol(species.name)
        self.name = name
        self.initial_count = 0
        self.color = color

    def __str__(self):
        return f'Marker with name {self.name} for species {repr(self.species)}'

    def __repr__(self):
        return f'{self.__class__.__name__}(name={repr(self.name)}, species={repr(self.species)}'