""""""
import math
from typing import Optional, Set, Mapping, Dict, Tuple, List

import logging
import copy

import sympy as sym
from monty.json import jsanitize, MontyDecoder
from sympy.parsing import sympy_parser

from dtcs.common import util
from dtcs.common.display import latex_map
from dtcs.spec.spec_abc import SpecCollection
from dtcs.spec.crn.sym_abc import SymSpec, ChemInfo


_logger = logging.getLogger(__name__)


class RxnABC(SymSpec):
    """A chemical reaction with reactants, products, and a rate constant.

    Attributes:
        reactants: A sym.Expr of the reactants.
        products: A sym.Expr of the products.
        rate_constant: A float, the rate constant of the chemical reaction.
    """

    def __init__(self, reactants: Optional[sym.Expr],
                 products: Optional[sym.Expr], k: float, **kwargs):
        """Create a new reaction by giving equation of the reactants.

        This is intended to look like reactants -> products @ rate k. That is,
        if your chemical equation looks like 2A -> B + C, then your reactants
        is 2A.

        Args:
            reactants: The left-hand side of the chemical reaction.
            products: The right-hand side of the chemical reaction.
            k: The rate constant.
        """
        super().__init__(**kwargs)
        self.rate_constant = k

        # Note that the type suggestion is Optional[sym.Expr].
        # It is possible that the user could pass None or 0 or 1 in.
        # Hence, sanitize input
        if not isinstance(reactants, sym.Expr):
            self.reactants = sym.sympify(0)
        else:
            self.reactants = reactants

        if not isinstance(products, sym.Expr):
            self.products = sym.sympify(0)
        else:
            self.products = products

    def get_symbols(self) -> Set[sym.Symbol]:
        symbols = set()
        symbols.update(self.reactants.free_symbols)
        symbols.update(self.products.free_symbols)
        return symbols

    @property
    def rate(self):  # TODO(Andrew): Remove, for bodge
        return self.rate_constant

    @property
    def name(self):
        return str(self)

    @property
    def reactant_symbols(self):
        """
        A list of free symbols.

        :return: a list of symbols on the reactant side.
        """
        return [s for s in self.reactants.free_symbols]

    @property
    def product_symbols(self):
        """
        A list of free symbols.

        :return: a list of symbols on the reactant side.
        """
        return [s for s in self.products.free_symbols]

    def rename(self, mapping: Mapping):
        self.reactants.subs(mapping)
        self.products.subs(mapping)

    def as_dict(self, sanitize=True) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict(sanitize=False)
        d['reactants'] = str(d['reactants'])
        d['products'] = str(d['products'])
        if sanitize:
            d = jsanitize(d, strict=True)
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        d['reactants'] = sympy_parser.parse_expr(d['reactants'])
        d['products'] = sympy_parser.parse_expr(d['products'])
        d['k'] = d.pop('rate_constant')
        return super(RxnABC, cls).from_dict(d)

    @util.depreciate
    def text(self) -> str:
        reactants: Dict[str, int] = self.reactants.as_coefficients_dict()
        products: Dict[str, int] = self.products.as_coefficients_dict()

        text = self.coeff_dict_to_text(reactants)

        conversion_verb = "are"
        if len(reactants) == 1:
            conversion_verb = "is"
        text += f" {conversion_verb} converted to "
        text += self.coeff_dict_to_text(products)

        text += f" at a rate of {self.rate_constant}."
        return text

    @util.depreciate
    def coeff_dict_to_text(self, coeff_dict) -> str:
        """Given a coefficient dictionary, convert it to a string.
        """
        items, text = coeff_dict.items(), ""
        for i, item in enumerate(items):
            symbol, coeff = item
            conjunction, units, end_comma = "", "units", ","
            if len(items) <= 2:
                end_comma = ""
            if i == len(items)-1 and len(items) > 1:
                conjunction, end_comma = "and ", ""
            if coeff == 1:
                units = "unit"
            text += f"{conjunction}{coeff} {units} of {symbol}{end_comma} "

        return text[:-1]

    def _repr_latex_(self) -> Optional[str]:
        return '$' + latex_map.sym_subs(self.reactants) \
               + r' \longrightarrow ' \
               + latex_map.sym_subs(self.products) \
               + r'\ \ @ k\! = \! ' + f'{self.rate_constant:.3f}' \
               + '$'

    def __str__(self):
        return f'{self.reactants} -> {self.products} @ k={self.rate_constant}'

    def __repr__(self):
        return f'{self.__class__.__name__}' \
               f'(reactants={repr(self.reactants)}, ' \
               f'products={repr(self.products)}, ' \
               f'k={self.rate_constant})'

    @util.depreciate
    def id(self):
        """Return a unique identifier for this reaction
        """
        return [f'{self.reactants}->{self.products}@{self.rate_constant}']

    @util.depreciate
    def fingerprint(self):
        """Return a unique identifier for this reaction, ignoring the reaction constant."""
        return [f'{self.reactants}->{self.products}']

    @util.depreciate
    def set_rate(self, rate: float) -> None:
        """
        Set or reset the rate constant in this reaction.

        :param rate: the new rate constant.
        :return: None
        """
        self.rate_constant = rate


class RevRxnABC(RxnABC):
    """A reversible reaction, essentially a reaction with two rate constants.

    Its use is to be quickly unpacked into two Rxns.
    """

    def __init__(self, reactants: Optional[sym.Expr],
                 products: Optional[sym.Expr], k: float,
                 k2: Optional[float] = None, **kwargs):
        """Create a reversible reaction by giving equation.

        This is intended to look like reactants <-> products @ rate k1,
        with the reverse rate k2.

        Args:
            reactants: The left-hand side of the chemical reaction.
            products: The right-hand side of the chemical reaction.
            k: The rate constant.
            k2: Optional, the rate constant for the reverse reaction. If not
                supplied, it's assumed to be 1/k.
        """

        super().__init__(reactants=reactants, products=products, k=k, **kwargs)
        self.rate_constant_reverse = k2 or 1 / self.rate_constant

    def to_rxns(self) -> Tuple[RxnABC, RxnABC]:
        return RxnABC(self.reactants, self.products, k=self.rate_constant), \
               RxnABC(self.products, self.reactants, k=self.rate_constant_reverse)

    @classmethod
    def from_dict(cls, d: dict):
        d['k2'] = d.pop('rate_constant_reverse')
        return super(RevRxnABC, cls).from_dict(d)

    def _repr_latex_(self) -> Optional[str]:
        latex = '$' + latex_map.sym_subs(self.reactants) \
               + r' \longleftrightarrow ' \
               + latex_map.sym_subs(self.products) \
               + r'\ \ @ k\! = \! ' + f'{self.rate_constant:.3f}' \

        if not math.isclose(
                self.rate_constant,
                1 / self.rate_constant_reverse,
                rel_tol=1e-3
        ):
            latex += r', k_-\! = \! ' + f'{self.rate_constant_reverse:.3f}' \

        latex += '$'
        return latex

    def __str__(self):
        return f'{self.reactants} <-> {self.products} ' \
               f'@ k={self.rate_constant}, k2={self.rate_constant_reverse}'

    def __repr__(self):
        return f'{self.__class__.__name__}' \
               f'(reactants={repr(self.reactants)}, ' \
               f'products={repr(self.products)}, ' \
               f'k={self.rate_constant}, k2={self.rate_constant_reverse})'

    def set_rates(self, rate: float, rate_reverse: Optional[float]) -> None:
        """
        Set the rate and its reverse.

        :param rate: the forward rate;
        :param rate_reverse: the reverse direction rate;
        :return: None.
        """
        self.set_rate(rate)
        if rate_reverse is None:
            rate_reverse = 1 / rate
        self.rate_constant_reverse = rate_reverse

    def id(self):
        """Return a unique identifier for this reaction
        This identifier is intentionally designed to return the same result as if
        same reversible reaction formed with two Rxns.
        """
        return [f'{self.reactants}->{self.products}@{self.rate_constant}',f'{self.products}->{self.reactants}@{self.rate_constant_reverse}']

    def fingerprint(self):
        """Return a unique identifier for this reaction, ignoring the reaction constants.

        This fingerprint is intentionally designed to return the same result as fingerprinting the
        same reversible reaction formed with two Rxns.
        """
        return [f'{self.reactants}->{self.products}',f'{self.products}->{self.reactants}']


class RxnSystemABC(SymSpec, SpecCollection):
    """TODO"""

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
            if not isinstance(component, (RxnABC, ChemInfo)):
                raise TypeError(f'Invalid component (positional argument) '
                                f'type {component.__class__.__name__} for '
                                f'{self.__class__.__name__}.')
            self.elements.append(component)

        # Pick an order for the symbols
        symbols = set()
        for component in self.elements:
            symbols.update(component.get_symbols())
        self._species = sorted([str(symbol) for symbol in symbols])

        # TODO: Remove this when possible
        symbol_index = {symbol: index for index, symbol in
                        enumerate(sorted(list(symbols), key=lambda s: str(s)))}
        self._symbol_index: Dict[sym.Symbol, int] = symbol_index

    @property
    def reactions(self) -> List[RxnABC]:
        return self.by_subclass()[RxnABC]

    def get_symbols(self) -> Set[sym.Symbol]:
        return set(self.symbol_index.keys())

    @property
    @util.depreciate
    def symbol_index(self):
        return self._symbol_index

    @util.depreciate
    def get_symbols_ordered(self) -> List[sym.Symbol]:
        symbols = [None] * len(self.get_symbols())
        for symbol in self.get_symbols():
            symbols[self.symbol_index[symbol]] = symbol
        return symbols

    @property
    def species(self) -> List[str]:
        return self._species

    @property
    def species_symbols(self):
        return [sym.Symbol(name) for name in self.species]

    def get_rates(self):
        rates = []
        for reaction in self.reactions:
            if isinstance(reaction, RevRxnABC):
                rates.append((reaction.rate_constant, reaction.rate_constant_reverse))
            else:
                rates.append(reaction.rate_constant)
        return rates

    def subs_rates(self, rates):
        """
        TODO(Andrew)
        """
        rsys = copy.deepcopy(self)
        index = 0
        for elem in rsys.elements:
            if isinstance(elem, RevRxnABC):
                if isinstance(rates[index], tuple):
                    elem.rate_constant = rates[index][0]
                    elem.rate_constant_reverse = rates[index][1]
                else:
                    elem.rate_constant = rates[index]
                    elem.rate_constant_reverse = 1 / rates[index]
                index += 1
            elif isinstance(elem, RxnABC):
                elem.rate_constant = rates[index]
                index += 1
        return rsys

    def rename(self, mapping: Mapping):
        for comp in self.elements:
            comp.rename(mapping)

        symbol_index = {}
        for symbol, index in self.symbol_index.items():
            symbol_index[symbol.subs(mapping)] = index
        self.symbol_index = symbol_index

    def insert(self, index, component):

        if not isinstance(component, (RxnABC, ChemInfo)):
            raise TypeError(f'Invalid component (positional argument) '
                            f'type {component.__class__.__name__} for '
                            f'{self.__class__.__name__}.')
        if index != len(self.elements):
            raise NotImplementedError('RxnSystem currently does not support'
                                      'insertion.')
        self.elements.append(component)

        # Remake self._species and the symbol index.
        _logger.warning('Remaking _species and _symbol_index, this is unstable '
                        'code.')
        symbols = set()
        for component in self.elements:
            symbols.update(component.get_symbols())
        self._species = sorted([str(symbol) for symbol in symbols])

        # TODO: Remove this when possible
        symbol_index = {symbol: index for index, symbol in
                        enumerate(sorted(list(symbols), key=lambda s: str(s)))}
        self._symbol_index: Dict[sym.Symbol, int] = symbol_index

    @classmethod
    def from_dict(cls, d: dict):
        decode = MontyDecoder().process_decoded
        d['components'] = [decode(component) for component in d['components']]
        d['symbol_index'] = {sym.Symbol(name): index for name, index in
                             d['symbol_index'].items()}
        return super(RxnSystemABC, cls).from_dict(d)

    def _repr_latex_(self) -> Optional[str]:
        latex = r'$ \text{Reaction System: } \newline \begin{gathered}'
        for index, rxn in enumerate(self.reactions):
            latex += f'({index + 1}) & ' \
                     + rxn._repr_latex_()[1:-1] \
                     + r' \newline '
        latex += r' \end{gathered} $'
        return latex

    def __str__(self):
        s = self.__class__.__name__ + ' with components:\n'
        for component in self.elements:
            comp_lines = str(component).splitlines()
            s += ''.join([f'\t{line}\n' for line in comp_lines])
        return s[:-1]

    def __repr__(self):
        return f'{self.__class__.__name__}(components={repr(self.elements)})'
