"""TODO(Andrew)"""

from typing import Mapping, Set

import copy

import sympy as sym
from sympy.parsing import sympy_parser
from monty.json import jsanitize

from lblcrn.spec.crn.sym_abc import ChemInfo


class Coverage(ChemInfo):
    """TODO(Andrew)"""

    def __init__(self,
                 site: sym.Symbol,
                 species: sym.Symbol,
                 coverage: float,
                 **kwargs):
        super().__init__(symbol=species, **kwargs)

        # TODO(Andrew) This is stupid, each species has exactly one site,
        #  otherwise it's a different species!
        self._site = site
        self.coverage = coverage
        if not (0 <= coverage <= 1):
            raise ValueError('Coverage should be between 0 and 1.')

    def get_symbols(self) -> Set[sym.Symbol]:
        return {self.symbol, self.site}

    @property
    def species(self):
        return str(self.symbol)

    @property
    def site(self):
        # TODO(Andrew) Maybe we should have species and _species, to match
        #  this design pattern
        return str(self._site)

    def rename(self, mapping: Mapping):
        self.symbol.subs(mapping)
        self._site.subs(mapping)

    def as_dict(self, sanitize=True) -> dict:
        """Return a MSON-serializable dict representation."""
        d = super().as_dict(sanitize=False)
        d['site'] = str(d['_site'])
        d['species'] = str(d['species'])
        if sanitize:
            d = jsanitize(d)
        return d

    @classmethod
    def from_dict(cls, d: dict):
        """Load from a dict representation."""
        # TODO(Andrew) Maybe have the __init__ accept string and symbol, so
        #  that sympy is optional? Then it'd be easier to serialize too.
        d['site'] = sympy_parser.parse_expr(d['site'])
        d['species'] = sympy_parser.parse_expr(d['species'])
        return super(ChemInfo, cls).from_dict(d)
