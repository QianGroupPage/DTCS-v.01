"""TODO"""

from typing import Dict, List, Mapping, Set
import sympy as sym

from monty.json import MontyDecoder

from lblcrn.common import util
from lblcrn.spec.crn.crn_abc import SymSpec, ChemInfo
from lblcrn.spec.crn.reaction import Rxn


class RxnSystem(SymSpec):
    """TODO"""

    _schema = [
        'symbol_index'
    ]

    _default = {
        'components': list,
    }

    def __init__(self, *components, **kwargs):
        """TODO: Old

        Create a new reaction system. Requires a SpeciesManager.

        Accepts Rxns, Revrxns, Concs, Schedules, Terms, ConcEqs,
        ConcDiffEqs, and (one) SpeciesManager in any order.

        If you have a function returning a collection of the above, you do
        not have to worry about unpacking the collection: it will unpack and
        flatten lists and tuples for you.
        """
        super().__init__(**kwargs)
        # Flatten so that you can give lists/tuples of components
        components = list(util.flat(components))

        for component in components:
            # Don't accept bad types
            if not isinstance(component, (Rxn, ChemInfo)):
                raise TypeError(f'Invalid component (positional argument) '
                                f'type {component.__class__.__name__} for '
                                f'{self.__class__.__name__}.')
            self.components.append(component)

        # Pick an order for the symbols
        symbols = set()
        for component in self.components:
            symbols.update(component.get_symbols())
        symbol_index = {symbol: index for index, symbol in
                        enumerate(sorted(list(symbols), key=lambda s: str(s)))}
        self.symbol_index: Dict[sym.Symbol, int] = symbol_index

    def get_symbols(self) -> Set[sym.Symbol]:
        return set(self.symbol_index.keys())

    def get_symbols_ordered(self) -> List[sym.Symbol]:
        symbols = [None] * len(self.get_symbols())
        for symbol in self.get_symbols():
            symbols[self.symbol_index[symbol]] = symbol
        return symbols

    def rename(self, mapping: Mapping):
        for comp in self.components:
            comp.rename(mapping)

        symbol_index = {}
        for symbol, index in self.symbol_index.items():
            symbol_index[symbol.subs(mapping)] = index
        self.symbol_index = symbol_index

    @classmethod
    def from_dict(cls, d: dict):
        decode = MontyDecoder().process_decoded
        d['components'] = [decode(component) for component in d['components']]
        d['symbol_index'] = {sym.Symbol(name): index for name, index in
                             d['symbol_index'].items()}
        return super(RxnSystem, cls).from_dict(d)

    def __str__(self):
        s = self.__class__.__name__ + ' with components:\n'
        for component in self.components:
            comp_lines = str(component).splitlines()
            s += ''.join([f'\t{line}\n' for line in comp_lines])
        return s[:-1]

    def __repr__(self):
        return f'{self.__class__.__name__}(components={repr(self.components)})'

