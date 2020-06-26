from lblcrn.crn_sym import Rxn, RevRxn, Surface, SurfaceRxn, SurfaceRevRxn, Schedule, Conc
import re
import random


def to_str(rule):
    """
    :param rule: a reaction rule
    :return: a list of strings representing the rule, suitable for Surface CRN
    """
    if isinstance(rule, SurfaceRevRxn):
        left = rule.reactants_str
        right = rule.products_str
        return [f"({rule.rate_constant}) {left} -> {right}",
                f"({rule.rate_constant_reverse}) {right} -> {left}"]
    elif isinstance(rule, SurfaceRxn):
        return [f"({rule.rate_constant}) {rule.reactants_str} -> {rule.products_str}"]
    elif isinstance(rule, RevRxn):
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
    for c in rsys.schedules:
        # TODO
        if isinstance(c, Schedule) and not isinstance(c, Conc):
            raise Exception(f"{c} is a schedule, not an initial concentration." +
                            f" This is not currently supported in the Surface CRN.")

        # if not re.match('[1-9][0-9]*', str(c.expression)):
        #     raise Exception(f"{c.symbol} must have a positive integer number of initial counts. Currently, " +
        #                     f"it's initial count is {c.expression}")
        species.extend([str(c.symbol) for _ in range(int(str(c.concentration)))])

    if not rsys.surface:
        rsys.surface = Surface('surface', (10, 10), color=(34, 139, 34))
        print("Using a default surface")

    rows = rsys.surface.size[0]
    cols = rsys.surface.size[1]

    random.seed(random_seed)
    random_indices = random.sample(range(rows * cols), len(species))

    surface_strs = [[rsys.surface.name for _ in range(cols)] for _ in range(rows)]
    for i, choice in enumerate(random_indices):
        surface_strs[choice // cols][choice % cols] = species[i]

    return "\n".join([" ".join(strs) for strs in surface_strs]) + "\n"
