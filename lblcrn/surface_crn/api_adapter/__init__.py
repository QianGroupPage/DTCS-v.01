"""A module for adapting the reaction systems generated in crn_sym into
api used by the core surface_crn engine used in surface_crns
"""

from lblcrn.crn_sym import Rxn, RevRxn
import re


def to_str(rule):
    """
    :param rule: a reaction rule
    :return: a list of strings representing the rule, suitable for Surface CRN
    """
    
    if isinstance(rule, RevRxn):
        return [f"({rule.rate_constant}) {rule.reactants} -> {rule.products}",
                f"({rule.rate_constant_reverse}) {rule.reactants} -> {rule.products}"]
    elif isinstance(rule, Rxn):
        return [f"({rule.rate_constant}) {rule.reactants} -> {rule.products}"]
    else:
        raise Exception("Rule is not of type crn_sym.Rxn or crn_sym.RevRxn.")


def generate_rules(rsys):
    """
    :param rsys: a reaction system object
    :return: a list of strings representing the rule, suitable for Surface CRN    
    """
    reactions = []
    for c in rsys.components:
        if isinstance(c, Rxn):
            reactions.extend(to_str(c))
    return reactions


def generate_initial_surface(rsys, random_seed=30):
    """
    :rsys: a reaction system object
    :return: a string representing the initial surface.
    """
    species = []
    for c in rsys.conc_eqs:
        # TODO
        if not re.match('[1-9][0-9]*', str(c.expression)):
            raise Exception(f"{c.symbol} must have a positive integer number of initial counts. Currently, it's initial count is {}")
        species.extend([str(c.symbol) for _ in range(int(str(c.expression))])
            

    if not rsys.surface:
        rsys.surface = SurfaceManager('surface', rows=10, cols=10, color=(34,139,34))
        print("Using a default surface")

    rows = rsys.surface.rows
    cols = rsys.surface.cols
    
    random.seed(random_seed)
    random_indices = random.sample(range(rows * cols), len(species))
    
    surface_strs = [[rsys.surface.name for _ in range(cols)] for _ in range(rows)]
    for i, choice in enumerate(random_choices):
        surface_strs[choice // cols][choice % cols] = species[i] 
        
    return "\n".join([" ".join(strs) for strs in surface_strs]) + "\n"
