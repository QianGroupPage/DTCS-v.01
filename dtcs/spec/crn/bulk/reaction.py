"""Classes for defining a chemical reaction system.

Exports:
    Rxn: A chemical reaction
    RevRxn: A reversible chemical reaction
    RxnSystem: A collection of reactions and conditions (e.g. initial
        concentrations).

Usage:
    Usage essentially always looks like the following, where you make a
    reaction system with other classes as inputs.

    RxnSystem(
        sm,

        Rxn(x, y, k=3.2),
        RevRxn(x + y, 2z, k1=0.01, k2=100),

        Conc(x, 2),
        ...
    )

    It looks like this because it was made to imitate the predecessor
    mathematica project.
"""

from typing import List, Tuple

import sympy as sym

from dtcs.spec.crn.rxn_abc import RxnABC, RevRxnABC
from dtcs.spec.crn.bulk.conditions import Term


class BulkRxn(RxnABC):

    def to_terms_nok(self, k_name: str):
        # Get the lefts and right
        # Here we're assuming that reactants and products are linear.
        lefts = self.reactants.as_coefficients_dict()
        rights = self.products.as_coefficients_dict()

        # Make a dict (symbol : term), initializing each to 0.
        term_dict = {}
        for symbol in self.get_symbols():
            term_dict[symbol] = sym.sympify(0)

        # Make the common part of all the terms
        common_part = sym.Symbol(k_name)
        for symbol in self.reactants.free_symbols:
            common_part *= (symbol ** lefts[symbol])

        # Populate the lefts
        for symbol in self.reactants.free_symbols:
            term_dict[symbol] += -1 * lefts[symbol] * common_part

        # Populate the rights
        for symbol in self.products.free_symbols:
            term_dict[symbol] += rights[symbol] * common_part

        # Make it into terms
        terms = []
        for symbol in term_dict.keys():
            terms.append(Term(symbol, term_dict[symbol]))

        return terms

    # TODO(Andrew): To reduce duplicate code, maybe make this depend on
    #  to_terms_nok. There's just an issue with name conflicts.
    def to_terms(self) -> List[Term]:
        """Create a list of terms from the reaction.

        Each term is essentially the reaction rate but positive or negative,
        depending on if it is a reactant or the product.

        This uses the reaction rate formula r =
        k*(concentration of products^coefficient of product):
        if this is false for your chemical equation, you will get an error.
        """

        # Get the lefts and right
        # Here we're assuming that reactants and products are linear.
        lefts = self.reactants.as_coefficients_dict()
        rights = self.products.as_coefficients_dict()

        # Make a dict (symbol : term), initializing each to 0.
        term_dict = {}
        for symbol in self.get_symbols():
            term_dict[symbol] = sym.sympify(0)

        # Make the common part of all the terms
        common_part = self.rate_constant
        for symbol in self.reactants.free_symbols:
            common_part *= (symbol ** lefts[symbol])

        # Populate the lefts
        for symbol in self.reactants.free_symbols:
            term_dict[symbol] += -1 * lefts[symbol] * common_part

        # Populate the rights
        for symbol in self.products.free_symbols:
            term_dict[symbol] += rights[symbol] * common_part

        # Make it into terms
        terms = []
        for symbol in term_dict.keys():
            terms.append(Term(symbol, term_dict[symbol]))

        return terms

class BulkRevRxn(BulkRxn, RevRxnABC):

    def to_rxns(self) -> Tuple[BulkRxn, BulkRxn]:
        return BulkRxn(self.reactants, self.products,
                       k=self.rate_constant), \
               BulkRxn(self.products, self.reactants,
                       k=self.rate_constant_reverse)

    def to_terms_nok(self, k_name: str):
        rxns = self.to_rxns()
        return [*rxns[0].to_terms_nok(k_name),
                *rxns[1].to_terms_nok(f'{k_name}^*')]

    def to_terms(self) -> List[Term]:
        rxns = self.to_rxns()
        return [*rxns[0].to_terms(), *rxns[1].to_terms()]
