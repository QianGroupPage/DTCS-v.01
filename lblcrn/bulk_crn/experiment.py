"""TODO

"""

from typing import List, Optional, Union

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
    def species(self) -> List[sym.Symbol]:
        """All the species in the experiment."""
        return []

    # --- Plotting -----------------------------------------------------------

    def plot(self, ax: Optional[plt.Axes] = None,
             species: Union[sym.Symbol, List[sym.Symbol], None] = None,
             ignore: Union[sym.Symbol, List[sym.Symbol], None] = None,
             **kwargs):
        """Plot the experiment.

        Args:
            ax: Optional, the axis on which to plot.
            species: Optional, the species you want to plot. Defaults to all.
            ignore: Optional, the species you don't want to plot. Defaults to
                none.
            **kwargs: Extra arguments, forwarded many places.
        """
        species = self._get_species_not_ignored(species, ignore)
        if ax is None:
            ax = plt.gca()

        self._plot(species=species, ax=ax)

    @abc.abstractmethod
    def _plot(self, ax: plt.Axes, species: List[sym.Symbol], **kwargs):
        """Do the actual plotting legwork.

        In general, this shouldn't call plt.show(), as that screws up axes.

        Args:
            ax: The plt.Axes on which to plot.
            species: A list of sym.Symbols, the species to plot.
            **kwargs: Forwarded.
        """
        pass

    # --- Utility ------------------------------------------------------------

    def _get_species_not_ignored(self,
                                 species: Union[sym.Symbol,
                                                List[sym.Symbol], None] = None,
                                 ignore: Union[sym.Symbol,
                                               List[sym.Symbol], None] = None):
        """Get the specified subset of species.

        Gets [species...] less [ignored...]. If speices is None, defaults to
        all species.

        Args:
            species: A list of (or single) symbol representing a species.
            ignore: A list of (or single) symbol representing a species.

        Returns:
            Species - ignored, with species defaulting to self.species.
        """
        # Wrap everything to be a list
        if isinstance(species, sym.Symbol):
            species = [species]
        if isinstance(species, sym.Symbol):
            species = [species]

        # Default values
        if not species:
            species = self.species
        if not ignore:
            ignore = []

        return [specie for specie in species if specie not in ignore]

    def _repr_html_(self) -> Optional[str]:
        """For iPython; mostly just gives the DataFrame."""
        # TODO: do this without accessing private member?
        df_html = self.df.head()._repr_html_()

        html = f"""<pre>{self.__class__.__name__} with head:</pre>
        {df_html}"""

        return html
