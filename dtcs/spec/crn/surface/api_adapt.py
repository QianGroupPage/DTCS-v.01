import random

from dtcs.sim.surface_crn.hex_grid_with_intersect import HexGridPlusIntersections
from dtcs.sim.surface_crn.surface_crns.options.option_processor import SurfaceCRNOptionParser
from dtcs.sim.surface_crn.surface_crns.readers.manifest_readers import read_manifest
from dtcs.spec.crn.bulk import Conc, RevRxn, Rxn, Schedule
from dtcs.spec.crn.surface.reaction import SurfaceRevRxn, SurfaceRxn

from dtcs.common.colors import color_map


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


def generate_initial_surface(surface, random_seed=30):
    """
    Generate an initial text-based surface.

    :rsys: a rxn_system object
    :return: a string representing the initial surface.
    """
    rows = surface.size[0]
    cols = surface.size[1]

    surface_strs = [[surface.name for _ in range(cols)] for _ in range(rows)]
    return 'UNUSED (Andrew)' + "\n".join([" ".join(strs) for strs in surface_strs]) + "\n"

    # species = []
    # locs = {}
    # for c in rsys.schedules:
    #     if isinstance(c, Schedule) and not isinstance(c, Conc):
    #         raise Exception(f"{c} is a schedule, not an initial concentration." +
    #                         f" This is not currently supported in the Surface CRN.")
    #
    #     if c.locs:
    #         locs[str(c.symbol)] = c.locs
    #     num_random = int(str(c.concentration)) - len(c.locs)
    #
    #     # if not re.match('[1-9][0-9]*', str(c.expression)):
    #     #     raise Exception(f"{c.symbol} must have a positive integer number of initial counts. Currently, " +
    #     #                     f"it's initial count is {c.expression}")
    #     species.extend([str(c.symbol) for _ in range(int(str(c.concentration)))])
    #
    # if not rsys.surface:
    #     rsys.surface = Surface('surface', (10, 10), color=(34, 139, 34))
    #     print("Using a default surface")
    #
    # rows = rsys.surface.size[0]
    # cols = rsys.surface.size[1]
    #
    # surface_strs = [[rsys.surface.name for _ in range(cols)] for _ in range(rows)]
    # if rsys.surface.use_coord_grid:
    #     surface_comments = "# The following surface is an unused dummy; a CoordGrid has been used."
    # elif rsys.surface.structure == "hexagon":
    #     surface_comments = "# The following surface is an unused dummy; a hexagonal grid has been used."
    # else:
    #     surface_comments = ""
    #
    #     random.seed(random_seed)
    #     random_indices = random.sample(range(rows * cols), len(species))
    #
    #     for i, choice in enumerate(random_indices):
    #         surface_strs[choice // cols][choice % cols] = species[i]
    # return surface_comments + "\n".join([" ".join(strs) for strs in surface_strs]) + "\n"

def generate_surface(rsys, random_seed=30):
    if rsys.surface is None:
        return None

    if rsys.surface.use_coord_grid:
        #  TODO: this is very horrible design!
        rsys.surface.set_default_name("top", rsys.surface.name)
        initial_concentrations = {site_name: {} for site_name in rsys.surface.allowed_sites}

        for c in rsys.schedules:
            # TODO
            if isinstance(c, Schedule) and not isinstance(c, Conc):
                raise Exception(f"{c} is a schedule, not an initial concentration." +
                                f" This is not currently supported in the Surface CRN.")

            if not rsys.species_manager.species_from_symbol(c.symbol).site or \
                   rsys.species_manager.species_from_symbol(c.symbol).site.name == "top":
                initial_concentrations["top"][str(c.symbol)] = int(str(c.concentration))
            else:
                initial_concentrations[rsys.species_manager.species_from_symbol(c.symbol).site.name][str(c.symbol)] = \
                        int(str(c.concentration))

            rsys.surface.set_initial_concentrations(initial_concentrations)

            # TODO: rename the sites in coord_grid
        return None
    elif rsys.surface.structure == "hexagon":
        species = {"threefold": [], "top": []}
        for c in rsys.schedules:
            # TODO
            if isinstance(c, Schedule) and not isinstance(c, Conc):
                raise Exception(f"{c} is a schedule, not an initial concentration." +
                                f" This is not currently supported in the Surface CRN.")

            if not rsys.species_manager.species_from_symbol(c.symbol).site or \
                   rsys.species_manager.species_from_symbol(c.symbol).site.name == "top":
                species["top"].extend([str(c.symbol) for _ in range(int(str(c.concentration)))])
            elif rsys.species_manager.species_from_symbol(c.symbol).site.name == "threefold":
                species["threefold"].extend([str(c.symbol) for _ in range(int(str(c.concentration)))])

        # Initiate the hexagonal grid
        rows = rsys.surface.size[0]
        cols = rsys.surface.size[1]
        surface = HexGridPlusIntersections(rows, cols)
        for n in surface:
            if n.is_intersection:
                n.state = "threefold"
            else:
                n.state = rsys.surface.name

        # Popultate the structure with initial atoms.
        initial_nodes = {"threefold": [], "top": []}
        for n in surface:
            if n.is_intersection:
                initial_nodes["threefold"].append(n)
            else:
                initial_nodes["top"].append(n)

        random.seed(random_seed)
        for k, v in species.items():
            nodes = initial_nodes[k]
            indices = random.sample(range(len(nodes)), len(v))
            chosen = [nodes[i] for i in indices]
            for i, n in enumerate(chosen):
                n.state = species[k][i]

        return surface
    elif rsys.surface.structure == "rectangle":
        return None


def generate_settings(rsys, max_duration, random_seed=923123122, video_path="Surface CRN Videos"):
    """
    :rsys: a rxn_system object
    :return: a string representing the initial surface.
    """
    return f"# Run settings\n" + \
           "pixels_per_node     = 80\n" + \
           "speedup_factor      = 0.5\n" + \
           "debug               = False\n" + \
           f"rng_seed           = {random_seed}\n" + \
           f"max_duration       = {max_duration}\n" + \
           "fps                 = 1\n" + \
           "node_display        = text\n" + \
           "wrap                = false\n" + \
           f"capture_directory  = {video_path}\n" + \
           "movie_title = SCRN Simulation\n\n"


def generate_colors(rsys):
    color_strs = ""
    for s in rsys.species:
        color = color_map.rgb256(s)
        color_strs += s + ": " + str(color) + "\n"
    return f"""!START_COLORMAP\n{color_strs}!END_COLORMAP\n"""


# TODO
def generate_manifest_stream(rsys, max_duration, random_seed_scrn=923123122, random_seed_surface=30,
                             video_path=""):
    """
    :param rsys: the rxn_system object
    :return: a stream of lines for the corresponding reaction rules.
    """
    rule = generate_settings(rsys, max_duration, random_seed_scrn, video_path=video_path)

    rule += "!START_TRANSITION_RULES\n"
    rule += "\n".join(generate_rules(rsys)) + "\n"
    rule += "!END_TRANSITION_RULES\n"
    rule += "\n"

    rule += "!START_INIT_STATE\n"
    rule += generate_initial_surface(rsys.surface, random_seed_surface)
    rule += "!END_INIT_STATE\n"
    rule += "\n"

    rule += generate_colors(rsys)

    for line in rule.splitlines():
        yield line


def _get_opts_via_manifest(rsys, time=0, rng_seed=0):
    manifest = generate_manifest_stream(
        rsys=rsys,
        max_duration=time,
        random_seed_scrn=rng_seed,
        random_seed_surface=0,  # This isn't being used
    )

    manifest_opts = read_manifest(
        filename=manifest,
    )

    opts = SurfaceCRNOptionParser(
        options=manifest_opts
    )

    return opts  # .__dict__