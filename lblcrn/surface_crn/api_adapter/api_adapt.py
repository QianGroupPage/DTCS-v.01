from lblcrn.crn_sym import Rxn, RevRxn, Surface, SurfaceRxn, SurfaceRevRxn, Schedule, Conc
from lblcrn.common import color_to_RGB
import re
import io
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
    :rsys: a rxn_system object
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


def generate_settings(rsys, max_duration, random_seed=923123122):
    """
    :rsys: a rxn_system object
    :return: a string representing the initial surface.
    """
    return f"# Run settings\n" + \
           "pixels_per_node     = 50\n" + \
           "speedup_factor      = 0.5\n" + \
           "debug               = False\n" + \
           f"rng_seed           = {random_seed}\n" + \
           f"max_duration       = {max_duration}\n" + \
           "fps                 = 1\n" + \
           "node_display        = text\n" + \
           "wrap                = false\n" + \
           f"capture_directory  = Surface CRN Videos\n" + \
           "movie_title = SCRN Simulation\n\n"


def generate_colors(rsys):
    color_strs = ""
    for s, color in rsys.get_colors().items():
        if isinstance(color, str):
            color = (c for c in color_to_RGB(color))
        color_strs += str(s) + ": " + str(color) + "\n"

    return f"""!START_COLORMAP\n{color_strs}!END_COLORMAP\n"""


# TODO
def generate_manifest_stream(rsys, max_duration, random_seed_scrn=923123122, random_seed_surface=30):
    """
    :param rsys: the rxn_system object
    :return: a stream of lines for the corresponding reaction rules.
    """
    rule = generate_settings(rsys, max_duration, random_seed_scrn)

    rule += "!START_TRANSITION_RULES\n"
    rule += "\n".join(generate_rules(rsys)) + "\n"
    rule += "!END_TRANSITION_RULES\n"
    rule += "\n"

    rule += "!START_INIT_STATE\n"
    rule += generate_initial_surface(rsys, random_seed_surface)
    rule += "!END_INIT_STATE\n"
    rule += "\n"

    rule += generate_colors(rsys)

    for line in rule.splitlines():
        yield line
