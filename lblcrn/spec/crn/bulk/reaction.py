"""TODO"""

from typing import List, Tuple

import sympy as sym

from lblcrn.spec.crn.reaction import Rxn, RevRxn
from lblcrn.spec.crn.bulk.conditions import Term


class BulkRxn(Rxn):

    def to_terms(self) -> List[Term]:
        """Create a list of terms from the reaction.

        Each term is essentially the reaction rate but positive or negative,
        depending on if it is a reactand or the product.

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

class BulkRevRxn(BulkRxn, RevRxn):

    def to_rxns(self) -> Tuple[BulkRxn, BulkRxn]:
        return BulkRxn(self.reactants, self.products,
                       k=self.rate_constant), \
               BulkRxn(self.products, self.reactants,
                       k=self.rate_constant_reverse)

    def to_terms(self) -> List[Term]:
        rxns = self.to_rxns()
        return [*rxns[0].to_terms(), *rxns[1].to_terms()]