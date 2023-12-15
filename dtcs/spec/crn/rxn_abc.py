""""""
import functools
from typing import Callable, Dict, List, Mapping, Optional, Set, Tuple, Union, TypeAlias
from abc import ABC, abstractmethod
from numbers import Number

import copy
import logging

import sympy as sym
from monty.json import jsanitize, MontyDecoder
from sympy.parsing import sympy_parser
from sympy.physics import units

from dtcs import config
from dtcs.common import util, const, display
from dtcs.spec.spec_abc import SpecCollection
from dtcs.spec.crn.sym_abc import SymSpec, ChemInfo, ChemExpression
from dtcs.common.const import RESERVED_SYMBOLS, DG, K, P, GIBBS_ENERGY, PRESSURE, TEMPERATURE

_logger = logging.getLogger(__name__)

Relation: TypeAlias = Union[sym.Expr, Callable, str, float]

DG_TO_RATE = {
    'basic': sym.exp(-0.1 * DG / (units.boltzmann * K)),
    'adsorption': P / (0.1 * units.torr) * \
                  sym.exp(-0.1 * DG / (units.boltzmann * K)),
}


class RxnABC(SymSpec, ABC):
    """A chemical reaction with reactants, products, and a rate constant.

    Attributes:
        reactants: A sym.Expr of the reactants.
        products: A sym.Expr of the products.
        rate_constant: A float, the rate constant of the chemical reaction.
    """

    def __init__(
            self,
            reactants: Optional[sym.Expr],
            products: Optional[sym.Expr], *,
            k: Optional[Relation] = None,
            dg: Optional[Relation] = None,
            **kwargs):

        """Create a new reaction by giving equation of the reactants.

        This is intended to look like reactants -> products @ rate k. That is,
        if your chemical equation looks like 2A -> B + C, then your reactants
        is 2A.

        Args:
            reactants: The left-hand side of the chemical reaction.
            products: The right-hand side of the chemical reaction.
            k: The rate constant.
            dg: The Gibbs free energy of the reaction.
        """
        super().__init__(**kwargs)

        # --- Class Variables ---
        self._gibbs_info: Optional[Relation]
        self._rate_info: Optional[Relation]
        self.reactants: sym.Expr
        self.products: sym.Expr

        # --- Initialize ---
        # Note that the type suggestion is Optional[sym.Expr].
        # It is possible that the user could pass None or 0 or 1 in.
        # Hence, sanitize input
        self.reactants = reactants if isinstance(reactants, sym.Expr) \
            else sym.sympify(0)
        self.products = products if isinstance(products, sym.Expr) \
            else sym.sympify(0)

        # If the rate is going to be calculated, we need gibbs energy
        if isinstance(k, (Callable, sym.Expr, str)) and (dg is None):
            raise TypeError('Cannot calculate rate without gibbs energy.')
        self._rate_info = k if (dg is None or k is not None) else 'basic'
        self._gibbs_info = dg

    # --- Properties ----------------------------------------------------------
    @property
    def name(self):  # TODO: Used internally, does it need to be unique?
        return str(self)

    @property
    def _is_fixed_rate(self) -> bool:
        return isinstance(self._rate_info, Number)

    @property
    def _is_func_rate(self) -> bool:
        return isinstance(self._rate_info, Callable)

    @property
    def _is_sym_rate(self) -> bool:
        return isinstance(self._rate_info, (sym.Expr, str))

    @property
    def _is_no_rate(self) -> bool:
        return self._rate_info is None

    # --- Chemistry -----------------------------------------------------------
    def get_rate(
            self,
            pressure: float = config.ref_tp['pressure'],
            temperature: float = config.ref_tp['temperature']
    ) -> Optional[float]:
        if self._is_fixed_rate:
            return float(self._rate_info)
        elif self._is_func_rate:
            return self._rate_info(self.get_gibbs(), pressure, temperature)
        elif self._is_sym_rate:
            dg_to_rate = _lambdify_dg_to_rate(self._rate_info)
            return dg_to_rate(self.get_gibbs(), pressure, temperature)

    def get_gibbs(self):
        return self._gibbs_info

    # --- Utility -----------------------------------------------------
    @abstractmethod
    def get_symbols(self) -> Set[sym.Symbol]:
        """Returns a set of all the sympy Symbols used in the reaction."""
        symbols = set()
        if self._is_sym_rate:
            symbols.update(_parse_sym_rate(self._rate_info).free_symbols)
        return symbols

    @abstractmethod
    def rename(self, mapping: Mapping):
        """According to the mapping, changes the names of all the sympy
        symbols or other relevant text-based components of the reaction."""
        if self._is_sym_rate:
            self._rate_info = _parse_sym_rate(self._rate_info).subs(mapping)

    # --- Serialization -------------------------------------------------------
    def as_dict(self, sanitize=True) -> dict:
        """Return a MSON-serializable dict representation."""
        raise NotImplementedError()
        d = super().as_dict(sanitize=False)
        d['reactants'] = str(d['reactants'])
        d['products'] = str(d['products'])
        if sanitize:
            d = jsanitize(d, strict=True)
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        raise NotImplementedError()
        d['reactants'] = sympy_parser.parse_expr(d['reactants'])
        d['products'] = sympy_parser.parse_expr(d['products'])
        d['k'] = d.pop('rate_constant')
        return super(RxnABC, cls).from_dict(d)

    # --- Representation ------------------------------------------------------
    @abstractmethod
    def _latex_rxn(self) -> str:
        """Get a LaTeX representation of the reaction, without worrying about
        the rate constants or gibbs energies. Don't add $ around it."""
        raise NotImplementedError()

    def latex(self, k_idx: Optional[int] = None) -> str:
        """Gets a LaTeX representation of the reaction, without $ around it."""
        k_forward, _ = const.k_names(k_idx)

        # Create a string for the gibbs energy beforehand
        unit_dg = sym.latex(config.units['energy'])
        dg_string = None if self.get_gibbs() is None else \
            r'\, @ \ \Delta G\! = \! ' + f'{self.get_gibbs():.3f}\\, {unit_dg}'

        # First we display the reaction
        st = self._latex_rxn()

        # Choose to display k or dG, whichever is a number
        if self._is_fixed_rate:
            st += rf'\, @ \ {k_forward}\! = \! {self.get_rate():.3f}'
        elif self._is_func_rate:
            st += dg_string + rf', \ {k_forward}\! \sim \! \Delta G'
        elif self._is_sym_rate:
            rate_info = _parse_sym_rate(self._rate_info).subs(const.PRETTY_SUBS)
            st += dg_string + rf', \ {k_forward}\! = \! {sym.latex(rate_info)}'

        return st

    def __str__(self):
        return f'{self.reactants} -> {self.products}'

    def __repr__(self):
        return f'{self.__class__.__name__}' \
               f'(reactants={repr(self.reactants)}, ' \
               f'products={repr(self.products)}, ' \
               f'k={repr(self._rate_info)}, ' \
               f'dg={repr(self._gibbs_info)})'

    # --- Depreciated ---------------------------------------------------------
    @property
    @util.depreciate
    def rate_constant(self):
        return self.get_rate()

    @property
    @util.depreciate
    def rate(self):  # TODO(Andrew): Remove, for bodge
        return self.rate_constant

    @property
    @util.depreciate
    def reactant_symbols(self):
        """
        A list of free symbols.

        :return: a list of symbols on the reactant side.
        """
        return [s for s in self.reactants.free_symbols]

    @property
    @util.depreciate
    def product_symbols(self):
        """
        A list of free symbols.

        :return: a list of symbols on the reactant side.
        """
        return [s for s in self.products.free_symbols]

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
            if i == len(items) - 1 and len(items) > 1:
                conjunction, end_comma = "and ", ""
            if coeff == 1:
                units = "unit"
            text += f"{conjunction}{coeff} {units} of {symbol}{end_comma} "

        return text[:-1]

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

    _rxn_cls = RxnABC

    def __init__(
            self,
            reactants: Optional[sym.Expr],
            products: Optional[sym.Expr],
            k: Optional[Relation] = None,
            k2: Optional[Relation] = None,
            dg: Optional[Relation] = None,
            **kwargs
    ):
        """Create a reversible reaction by giving equation.

        This is intended to look like reactants <-> products @ rate k,
        with the reverse rate k2.

        Args:
            reactants: The left-hand side of the chemical reaction.
            products: The right-hand side of the chemical reaction.
            k: The rate constant.
            k2: Optional, the rate constant for the reverse reaction. If not
                supplied, it's assumed to be 1/k.
            dg: The Gibbs free energy of the reaction.
        """

        super().__init__(
            reactants=reactants,
            products=products,
            k=k,
            dg=dg,
            **kwargs,
        )
        # --- Class Variables ---
        self._rev_rate_info: Optional[Relation]

        # --- Initialize ---
        if k2 is not None and self._is_no_rate:
            raise TypeError('Do not supply k2 without k.')
        elif isinstance(k2, (Callable, sym.Expr, str)) and (dg is None):
            raise TypeError('Cannot calculate rev. rate without gibbs energy.')

        self._rev_rate_info = k2

    # --- Properties ----------------------------------------------------------
    @property
    def _is_fixed_rev_rate(self) -> bool:
        return isinstance(self._rev_rate_info, Number)

    @property
    def _is_func_rev_rate(self) -> bool:
        return isinstance(self._rev_rate_info, Callable)

    @property
    def _is_sym_rev_rate(self) -> bool:
        return isinstance(self._rev_rate_info, (sym.Expr, str))

    @property
    def _is_no_rev_rate(self) -> bool:
        return self._rev_rate_info is None

    # --- Chemistry -----------------------------------------------------------
    def get_rev_rate(
            self,
            pressure: float = config.ref_tp['pressure'],
            temperature: float = config.ref_tp['temperature']
    ) -> Optional[float]:
        if self._is_fixed_rev_rate:
            return float(self._rev_rate_info)
        elif self._is_func_rev_rate:
            return self._rev_rate_info(-1 * self.get_gibbs(), pressure, temperature)
        elif self._is_sym_rev_rate:
            dg_to_rate = _lambdify_dg_to_rate(self._rev_rate_info)
            return dg_to_rate(-1 * self.get_gibbs(), pressure, temperature)
        # After we've tried using the supplied information about the reverse
        #  reaction, we default to trying to invert the forward reaction.
        elif self._is_fixed_rate:
            # If rate 1 exists we can default to 1 / k.
            return float(1 / self._rate_info)
        elif self._is_func_rate:
            return self._rate_info(-1 * self.get_gibbs(), pressure, temperature)
        elif self._is_sym_rate:
            dg_to_rate = _lambdify_dg_to_rate(self._rate_info)
            return dg_to_rate(-1 * self.get_gibbs(), pressure, temperature)

    def to_rxns(self) -> Tuple[_rxn_cls, _rxn_cls]:
        gibbs_info = self._gibbs_info and -1 * self._gibbs_info
        rev_rate_info = self._rev_rate_info if not self._is_no_rev_rate else (
            1 / self._rate_info if self._is_fixed_rate else self._rate_info
        )

        return self._rxn_cls(self.reactants, self.products,
                             k=self._rate_info, dg=self._gibbs_info), \
               self._rxn_cls(self.products, self.reactants,
                             k=rev_rate_info, dg=gibbs_info)

    # --- Utility -----------------------------------------------------
    def get_symbols(self) -> Set[sym.Symbol]:
        symbols = super().get_symbols()
        if self._is_sym_rev_rate:
            symbols.update(_parse_sym_rate(self._rev_rate_info).free_symbols)
        return symbols

    def rename(self, mapping: Mapping):
        super().rename(mapping)
        if self._is_sym_rev_rate:
            self._rev_rate_info = _parse_sym_rate(self._rev_rate_info).subs(mapping)

    # --- Serialization -----------------------------------------------------
    @classmethod
    def from_dict(cls, d: dict):
        raise NotImplementedError()
        d['k2'] = d.pop('rate_constant_reverse')
        return super(RevRxnABC, cls).from_dict(d)

    # --- Representation -----------------------------------------------------
    @abstractmethod
    def _latex_rxn(self) -> str:
        """Get a LaTeX representation of the reaction, without worrying about
        the rate constants or gibbs energies. Don't add $ around it."""
        raise NotImplementedError()

    def latex(self, k_idx: Optional[int] = None) -> str:
        """Gets a LaTeX representation of the reaction, without $ around it."""
        k_forward, k_reverse = const.k_names(k_idx)

        # Create a string for the gibbs energy beforehand
        unit_dg = sym.latex(config.units['energy'])
        dg_string = None if self.get_gibbs() is None else \
            r'\, @ \ \Delta G\! = \! ' + f'{self.get_gibbs():.3f}\\, {unit_dg}'

        # First we display the reaction
        st = self._latex_rxn()

        # Create the reverse rate's string in advance so we can add it by
        #  casework later and not have to worry we're doing the wrong one.
        rev_rate_str = ''
        if self._is_fixed_rev_rate:
            rev_rate_str = rf', \ {k_reverse}\! = \! {self.get_rev_rate():.3f}'
        elif self._is_sym_rev_rate:
            rev_rate_info = _parse_sym_rate(self._rev_rate_info).subs(const.PRETTY_SUBS)
            rev_rate_str = rf', \ {k_reverse}\! = \! {sym.latex(rev_rate_info)}'
        elif self._is_func_rev_rate:
            rev_rate_str = rf', \ {k_reverse}\! \sim \! \Delta G'

        if self._is_fixed_rate:
            st += rf'\, @ \ {k_forward}\! = \! {self.get_rate():.3f}'
            st += rev_rate_str
        elif self._is_func_rate and self._is_func_rev_rate:
            st += dg_string + rf', \ {k_forward}, {k_reverse}\! \sim \! \Delta G'
        elif self._is_func_rate:
            st += dg_string + rf', \ {k_forward}\! \sim \! \Delta G' + rev_rate_str
        elif self._is_sym_rate:
            rate_info = _parse_sym_rate(self._rate_info).subs(const.PRETTY_SUBS)
            st += dg_string
            st += rf', \ {k_forward}\! = \! {sym.latex(rate_info)}'
            st += rev_rate_str
        else:
            st += rev_rate_str

        return st

    def __str__(self):
        return f'{self.reactants} <-> {self.products}'

    def __repr__(self):
        return f'{self.__class__.__name__}' \
               f'(reactants={repr(self.reactants)}, ' \
               f'products={repr(self.products)}, ' \
               f'k={repr(self._rate_info)}, ' \
               f'k2={repr(self._rev_rate_info)}, ' \
               f'dg={repr(self._gibbs_info)})'

    # --- Depreciated ---------------------------------------------------------
    @property
    @util.depreciate
    def rate_constant_reverse(self):
        return self.get_rev_rate()

    @util.depreciate
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

    @util.depreciate
    def id(self):
        """Return a unique identifier for this reaction
        This identifier is intentionally designed to return the same result as if
        same reversible reaction formed with two Rxns.
        """
        return [f'{self.reactants}->{self.products}@{self.rate_constant}',
                f'{self.products}->{self.reactants}@{self.rate_constant_reverse}']

    @util.depreciate
    def fingerprint(self):
        """Return a unique identifier for this reaction, ignoring the reaction constants.

        This fingerprint is intentionally designed to return the same result as fingerprinting the
        same reversible reaction formed with two Rxns.
        """
        return [f'{self.reactants}->{self.products}', f'{self.products}->{self.reactants}']


class RxnSystemABC(SymSpec, SpecCollection):
    """TODO"""

    def __init__(self, *components, **kwargs):
        """Create a new reaction system.

        Accepts Rxns, Revrxns, Concs, Schedules, Terms, ConcEqs, ConcDiffEqs.

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
        # Remove non-species like time, pressure, and temperature.
        symbols -= RESERVED_SYMBOLS
        self._species = tuple(sorted([str(symbol) for symbol in symbols]))

    # --- Properties ----------------------------------------------------------
    @property
    def reactions(self) -> List[RxnABC]:
        # TODO: We are heavily trusting that this will not change.
        return self.by_subclass()[RxnABC]

    @property
    def network(self) -> List[Union[RxnABC, ChemExpression]]:
        return self.reactions + self.by_subclass()[ChemExpression]

    @property
    def species(self) -> Tuple[str]:
        # TODO: We are heavily trusting that this will not change.
        return self._species

    @property
    def rates(self) -> Tuple[str]:
        return tuple(self.get_rates().keys())

    @property
    def species_symbols(self) -> Tuple[sym.Symbol]:
        return tuple(sym.Symbol(name) for name in self.species)

    @property
    def rate_symbols(self) -> Tuple[sym.Symbol]:
        return tuple(sym.Symbol(name) for name in self.rates)

    # --- Chemistry -----------------------------------------------------------
    def get_rates(self,
                  pressure: float = config.ref_tp['pressure'],
                  temperature: float = config.ref_tp['temperature']):
        """Gets the rates from each of the reactions. This function relies on
        self.reactions not changing in order."""
        rates = {}
        for index, reaction in enumerate(self.reactions):
            k_forward, k_reverse = const.k_names(index)
            rates[k_forward] = reaction.get_rate(pressure, temperature)
            if isinstance(reaction, RevRxnABC):
                rates[k_reverse] = reaction.get_rev_rate(pressure, temperature)
        return rates

    def subs_gibbs(self, gibbs):
        """Returns a copy with new Gibbs energies from each of the reactions.
        This function relies on self.reactions not changing in order.
        """
        rsys = copy.deepcopy(self)
        for index, rxn in enumerate(rsys.reactions):
            rxn._gibbs_info = gibbs[index]
        return rsys

    def subs_rates(self, rates):
        """
        Returns a copy with new rates from each of the reactions. This function
        relies on self.reactions not changing in order.
        """
        rsys = copy.deepcopy(self)
        for index, rxn in enumerate(rsys.reactions):
            key_forward, key_reverse = const.k_names(index)

            # Simple part: set the forward rate
            rxn._rate_info = rates[key_forward]

            # Complicated part: set the reverse rate.
            #  We need to check if it's a reversible reaction, and if the
            #  user supplied rates[key_reverse] at all.
            #  Luckily, it's okay if ._rev_rate_info is None.
            if isinstance(rxn, RevRxnABC):
                rxn._rev_rate_info = rates.get(key_reverse)
        return rsys

    # --- Utility -------------------------------------------------------------
    def get_symbols(self) -> Set[sym.Symbol]:
        return set(self.symbol_index.keys())

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
        # Remove non-species like time, pressure, and temperature.
        symbols -= RESERVED_SYMBOLS
        self._species = sorted([str(symbol) for symbol in symbols])

    # --- Serialization -------------------------------------------------------
    @classmethod
    def from_dict(cls, d: dict):
        decode = MontyDecoder().process_decoded
        d['components'] = [decode(component) for component in d['components']]
        d['symbol_index'] = {sym.Symbol(name): index for name, index in
                             d['symbol_index'].items()}
        return super(RxnSystemABC, cls).from_dict(d)

    # --- Representation ------------------------------------------------------
    def latex(self) -> Optional[str]:
        latex = r'\text{Reaction System: } \newline \begin{gathered}'
        for index, elem in enumerate(self.network):
            # If it's a reaction, let's number the rates
            elem_latex = (
                elem.latex(k_idx=index)
                if isinstance(elem, RxnABC)
                else elem.latex()
            )
            latex += f'({index}) & ' \
                     + elem_latex \
                     + r' \newline '
        latex += r' \end{gathered}'
        return latex

    def __str__(self):
        s = self.__class__.__name__ + ' with components:\n'
        for component in self.elements:
            comp_lines = str(component).splitlines()
            s += ''.join([f'\t{line}\n' for line in comp_lines])
        return s[:-1]

    def __repr__(self):
        return f'{self.__class__.__name__}(components={repr(self.elements)})'

    # --- Depreciated ---------------------------------------------------------
    @property
    @util.depreciate
    def symbol_index(self):
        return {sym.Symbol(name): index for index, name in enumerate(self._species)}

    @util.depreciate
    def get_symbols_ordered(self) -> List[sym.Symbol]:
        symbols = [None] * len(self.get_symbols())
        for symbol in self.get_symbols():
            symbols[self.symbol_index[symbol]] = symbol
        return symbols


def _parse_sym_rate(expr_or_key: Union[sym.Expr, str]):
    """Given a sympy rate expression or string corresponding to a default,
    give a sympy rate expression via looking up any defaults."""
    if isinstance(expr_or_key, str):
        return DG_TO_RATE[expr_or_key]
    elif isinstance(expr_or_key, sym.Expr):
        return expr_or_key
    else:
        raise TypeError(f'Invalid type {type(expr_or_key)}')


@functools.lru_cache(maxsize=127)
def _lambdify_dg_to_rate(expr: Union[sym.Expr, str]):
    """Turns a sympy.Expr expression of DG, P, and K into a function."""
    expr = _parse_sym_rate(expr)
    gibbs, time, pressure, temp = sym.symbols('gibbs time pressure temp')

    dg_to_rate_expr = expr.subs({
        GIBBS_ENERGY: gibbs * config.units['energy'],
        PRESSURE: pressure * config.units['pressure'],
        TEMPERATURE: temp * config.units['temperature'],
    })

    dg_to_rate_expr = units.convert_to(dg_to_rate_expr, config.unit_system)

    return sym.lambdify(
        args=(gibbs, pressure, temp),
        expr=dg_to_rate_expr,
    )
