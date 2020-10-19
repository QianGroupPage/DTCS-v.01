"""TODO

"""

import abc
from typing import List, Optional, Union

import sympy as sym
from matplotlib import pyplot as plt

import lblcrn
import monty.json


class Experiment(monty.json.MSONable, abc.ABC):
    """A base class for experiments.

    Attributes:
        df: The pd.DataFrame representing the experiment.
    """

    def __init__(self):
        """Initialize required instance variables."""
        self.df = None

    # --- Plotting -----------------------------------------------------------

    def plot(self, ax: Optional[plt.Axes] = None,
             species: Union[sym.Symbol, List[sym.Symbol], None] = None,
             ignore: Union[sym.Symbol, List[sym.Symbol], None] = None,
             **kwargs) -> plt.Axes:
        """Plot the experiment.

        Args:
            ax: Optional, the axis on which to plot.
            species: Optional, the species you want to plot. Defaults to all.
            ignore: Optional, the species you don't want to plot. Defaults to
                none.
            **kwargs: Extra arguments, forwarded many places.

        Returns:
            #TODO
        """
        species = _get_species_not_ignored(species, ignore, self.species)
        if ax is None:
            ax = plt.gca()

        return self._plot(species=species, ax=ax, **kwargs)

    @abc.abstractmethod
    def _plot(self, ax: plt.Axes, species: List[sym.Symbol],
              **kwargs) -> plt.Axes:
        """Do the actual plotting legwork.

        In general, this shouldn't call plt.show(), as that screws up axes.

        Args:
            ax: The plt.Axes on which to plot.
            species: A list of sym.Symbols, the species to plot.
            **kwargs: Forwarded.

        Returns:
            #TODO
        """
        pass

    # --- Utility ------------------------------------------------------------

    @property
    @abc.abstractmethod
    def species(self) -> List[sym.Symbol]:
        """All the species in the experiment."""
        pass

    def as_dict(self) -> dict:
        """Return a MSON-serializable dict representation."""
        d = {
            '@module': self.__class__.__module__,
            '@class': self.__class__.__name__,
            '@version': lblcrn.__version__,  # TODO: Better way to do this?
        }
        return d

    @classmethod
    @abc.abstractmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        return cls()

    def _repr_html_(self) -> Optional[str]:
        """For iPython; mostly just gives the head of the DataFrame."""
        # TODO: do this without accessing private member?
        df_html = self.df.head()._repr_html_()

        html = f"""<pre>{self.__class__.__name__} with head:</pre>
        {df_html}"""

        return html

    def __str__(self) -> str:
        """Mostly just the head of the DataFrame."""
        # TODO: do this without accessing private member?
        df_str = str(self.df.head())

        s = f'{self.__class__.__name__} with head:\n{df_str}'

        return s

def _get_species_not_ignored(species: Union[sym.Symbol,
                                            List[sym.Symbol], None] = None,
                             ignore: Union[sym.Symbol,
                                           List[sym.Symbol], None] = None,
                             all_species: Union[sym.Symbol,
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
    if isinstance(ignore, sym.Symbol):
        ignore = [ignore]

    # Default values
    if not species:
        species = all_species
    if not ignore:
        ignore = []

    return [specie for specie in species if specie not in ignore]
