from lblcrn.surface_crn.surface_crns import SurfaceCRNQueueSimulator

from lblcrn.surface_crn.surface_crns.simulators.queue_simulator import *
from lblcrn.surface_crn.surface_crns.readers.manifest_readers import read_manifest
from lblcrn.surface_crn.surface_crns.options.option_processor import SurfaceCRNOptionParser
from lblcrn.surface_crn.results import Results
from lblcrn.common import color_to_HEX
from lblcrn.surface_crn.api_adapter.api_adapt import generate_manifest_stream, generate_surface,\
    HexGridPlusIntersectionDisplay
import os
from shutil import rmtree


def scrn_simulate(rxns, time_max=100, lattice=None, display_class=None, video=False, spectra_in_video=True,
                  running_average=10, species_tracked=[], manifest_file=""):
    """
    :param rxns:
    :param time_max:
    :param lattice:
    :param display_class:
    :param video:
    :param concentrations:
    :param species_tracked: a list of sympy symbols or strings representing each species
    :param manifest_file:
    :return:
    """
    if not manifest_file:
        manifest = generate_manifest_stream(rxns, time_max)
    else:
        manifest = manifest_file

    if not species_tracked:
        species_tracked = list(rxns.get_symbols())
    if not lattice:
        surface = generate_surface(rsys=rxns)
    else:
        surface = lattice
    #  TODO: infer spectra's scale from here.
    times, concs = simulate_without_display(manifest, surface, [str(s) for s in species_tracked], rxns)
    if manifest_file:
        r = Results.from_concs_times(manifest_file, rxns, concs, times)
    else:
        r = Results.from_concs_times("".join(list(manifest)), rxns, concs, times)

    video_link = None
    if video:
        # Generate the file stream again after it's used.
        if not manifest_file:
            manifest = generate_manifest_stream(rxns, time_max)
        video_link = get_video_link(manifest)
        # Generate the file stream again after it's used.
        if not manifest_file:
            manifest = generate_manifest_stream(rxns, time_max)
        # TODO: check to not overwrite the video files.
        frames_link = get_frames_link(manifest)
        if os.path.isdir(frames_link):
            rmtree(frames_link)
        # Generate the file stream again after it's used.
        if not manifest_file:
            manifest = generate_manifest_stream(rxns, time_max)
        if not lattice:
            surface = generate_surface(rsys=rxns)
        else:
            surface = lattice

        # TODO: fix the issue that the video will fail if there is already file in
        # the frames folder.
        # TODO: progress bar for the video
        # TODO: add this as an argument spectra_max_conc=r.df_raw.max()
        simulate_with_display(manifest, surface, rxns=rxns, spectra_in_video=spectra_in_video,
                              running_average=running_average, spectra_max_conc=r.df_raw.to_numpy().max())
    r.video = video_link

    # TODO: warn the user if termination is early.
    return r


def get_opts(manifest):
    '''
    Process a stream into an options file.
    '''
    manifest_options = read_manifest(manifest)
    return SurfaceCRNOptionParser(manifest_options)


def get_video_link(manifest):
    """
    :param manifest: a stream or file name of manifest file
    :return: the link where the video would be stored.
    """
    opts = get_opts(manifest)
    return f"{opts.capture_directory}/{opts.movie_title}.mp4"


# TODO: cleanup
def get_frames_link(manifest):
    opts = get_opts(manifest)
    return f"{opts.capture_directory}/frames"


def simulate_with_display(manifest_file, lattice, rxns=None, spectra_in_video=True, running_average=10,
                          spectra_max_conc=-1):
    if rxns.surface.structure == "hexagon":
        display_class = HexGridPlusIntersectionDisplay
    else:
        display_class = None
    SurfaceCRNQueueSimulator.simulate_surface_crn(manifest_file, display_class, init_state=lattice, rxns=rxns,
                                                  spectra_in_video=spectra_in_video, running_average=running_average,
                                                  spectra_max_conc=spectra_max_conc)


def add_groups(surface, rsys):
    """
    Update the surface with groups in the species manager.

    :param surface: a surface structure
    :param rsys: a rxn system
    """
    sm = rsys.species_manager
    size_dict = {s.name: s.size for s in sm.large_species}
    seen = set()
    for s in surface:
        if s not in seen and s.state in size_dict:
            group_size = size_dict[s.state]
            # Build the group
            free_neighbors = []
            for t in s.neighbors:
                n = t[0]
                if n not in seen and n.state in rsys.surface_names:
                    free_neighbors.append(n)
                seen.add(n)
            seen.add(s)

            # Pick from free neighbors as part of the group
            group = random.sample(free_neighbors, group_size - 1) + [s]

            for n in group:
                n.group = group
                n.state = s.state


def simulate_without_display(manifest_file, lattice, species_tracked, rxns):
    '''
    Run until completion or max time, storing an array of species counts at
    the times of each reaction, along with an array of times. At the end,
    display a graph of species concentrations as they change.
    '''
    manifest_options = read_manifest(manifest_file)

    opts = SurfaceCRNOptionParser(manifest_options)
    if not lattice:
        # If no grid is made, use the inåitial grid
        lattice = opts.grid

    # print(lattice)
    # TODO: something similar for the videos
    add_groups(lattice, rxns)
    # print(lattice)

    simulator = QueueSimulator(surface=lattice,
                               transition_rules=opts.transition_rules,
                               seed=opts.rng_seed,
                               simulation_duration=opts.max_duration)

    times = [0]
    concs = dict()
    for species in species_tracked:
        concs[species] = [0]
    for node in lattice:
        if node.state in concs:
            concs[node.state][0] += 1
    while not simulator.done():
        next_rxn = simulator.process_next_reaction()
        if next_rxn is None:
            break
        times.append(next_rxn.time)
        for species in species_tracked:
            concs[species].append(concs[species][-1])
        # Very simple mechanism, one reaction at a time.
        for reactant in next_rxn.rule.inputs:
            if reactant in concs:
                concs[reactant][-1] -= 1
        for product in next_rxn.rule.outputs:
            if product in concs:
                concs[product][-1] += 1

    return times, concs
