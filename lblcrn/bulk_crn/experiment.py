"""TODO

"""

from typing import List, Optional

import abc
from matplotlib import pyplot as plt
import sympy as sym


class Experiment(abc.ABC):
    """TODO

    """

    def __init__(self):
        """TODO

        """
        self.df = None

    @property
    @abc.abstractmethod
    def species(self):
        """TODO

        """
        pass

    def plot(self, show: bool = True,
             species: List[sym.Symbol] = [],
             ignore: List[sym.Symbol] = [],
             **kwargs):
        """TODO
        """
        species = self._get_species_not_ignored(species, ignore)

        self._plot(species)

        if show:
            plt.show()

    @abc.abstractmethod
    def _plot(self, species):  # TODO: Add type hints
        """TODO

        This method should _not_ call plt.show().

        Args:
            species:

        Returns:

        """
        pass

    # --- Utility ------------------------------------------------------------

    def _get_species_not_ignored(self, species: List[sym.Symbol] = None,
                                 ignore: List[sym.Symbol] = None):
        """A helper method to return species - ignored. If species is none,
        considers all species."""
        if not species:
            species = self.species
        if not ignore:
            ignore = []
        return [specie for specie in species if specie not in ignore]

    def _repr_html_(self) -> Optional[str]:
        """For iPython, mostly just gives the DataFrame."""
        # TODO: do this without accessing private member?
        df_html = self.df.head()._repr_html_()

        html = f"""<pre>{self.__class__.__name__} with head:</pre>
        {df_html}"""

        return html