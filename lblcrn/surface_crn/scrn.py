import os
from shutil import rmtree
from IPython.display import clear_output
from lblcrn.common import ipython_visuals
from lblcrn.surface_crn.api_adapter.api_adapt import generate_manifest_stream, \
    generate_surface, HexGridPlusIntersectionDisplay
from lblcrn.surface_crn.ensemble import Ensemble
from lblcrn.surface_crn.results import Results
from lblcrn.surface_crn.surface_crns import SurfaceCRNQueueSimulator
from lblcrn.surface_crn.surface_crns.simulators.queue_simulator import *
from lblcrn.surface_crn.surface_crns.readers.manifest_readers import read_manifest
from lblcrn.surface_crn.surface_crns.options.option_processor import SurfaceCRNOptionParser
from lblcrn.surface_crn.surface_crns.models.coord_grid import CoordGrid
from lblcrn.surface_crn.surface_crns.views.coord_grid_display import CoordGridDisplay


def scrn_simulate_single_run(rxns, time_max=100, lattice=None, display_class=None, video=False, spectra_in_video=True,
                             spectra_average_duration=2, species_tracked=[], manifest_file="", rng_seed=923123122,
                             video_path=""):
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
    # Dangerous! rng_seed + 1 is used
    group_selection_seed = rng_seed + 1

    if not manifest_file:
        manifest = generate_manifest_stream(rxns, time_max, random_seed_scrn=rng_seed, video_path=video_path)
    else:
        manifest = manifest_file

    if not species_tracked:
        species_tracked = list(rxns.get_symbols())

    # Provide a grid structure for use in place of the grid structure in the rules manifest.
    if not lattice and rxns.surface.use_coord_grid:
        # For CoordGrid, generate_surface function updates the surface object in accordance with
        # initial concentration and default species on the sites specified by the rules file.
        generate_surface(rsys=rxns)
        surface = CoordGrid.from_poscar(rxns.surface.poscar_file,
                                        supercell_dimensions=rxns.surface.supercell_dimensions,
                                        ignore_threhold=rxns.surface.surface_depth)
        surface.set_initial_concentrations(rxns.surface.initial_species)
        surface.set_default_species(rxns.surface.default_names)
    elif not lattice:
        surface = generate_surface(rsys=rxns)
    else:
        surface = lattice
    #  TODO: infer spectra's scale from here.
    times, concs = simulate_without_display(manifest, surface, [str(s) for s in species_tracked], rxns,
                                            group_selection_seed)

    # print(concs)

    if manifest_file:
        r = Results.from_concs_times(manifest_file, rxns, concs, times)
    else:
        r = Results.from_concs_times("".join(list(manifest)), rxns, concs, times)

    video_link = None
    if video:
        # Generate the file stream again after it's used.
        if not manifest_file:
            manifest = generate_manifest_stream(rxns, time_max, random_seed_scrn=rng_seed, video_path=video_path)
        video_link = get_video_link(manifest)

        # Generate the file stream again after it's used.
        if not manifest_file:
            manifest = generate_manifest_stream(rxns, time_max, random_seed_scrn=rng_seed, video_path=video_path)
        # TODO: check to not overwrite the video files.
        frames_link = get_frames_link(manifest)
        if os.path.isdir(frames_link):
            rmtree(frames_link)
        # Generate the file stream again after it's used.
        if not manifest_file:
            manifest = generate_manifest_stream(rxns, time_max, random_seed_scrn=rng_seed, video_path=video_path)
        # Provide a grid structure for use in place of the grid structure in the rules manifest.
        if not lattice and rxns.surface.use_coord_grid:
            # For CoordGrid, generate_surface function updates the surface object in accordance with
            # initial concentration and default species on the sites specified by the rules file.
            generate_surface(rsys=rxns)
            surface = CoordGrid.from_poscar(rxns.surface.poscar_file,
                                            supercell_dimensions=rxns.surface.supercell_dimensions,
                                            ignore_threhold=rxns.surface.surface_depth)
            surface.set_initial_concentrations(rxns.surface.initial_species)
            surface.set_default_species(rxns.surface.default_names)
        elif not lattice:
            surface = generate_surface(rsys=rxns)
        else:
            surface = lattice

        # TODO: fix the issue that the video will fail if there is already file in
        # the frames folder.
        # TODO: progress bar for the video
        # TODO: add this as an argument spectra_max_conc=r.df_raw.max()
        r.video_trajectory = simulate_with_display(manifest, surface, group_selection_seed, rxns=rxns,
                                                   spectra_in_video=spectra_in_video,
                                                   running_average=spectra_average_duration,
                                                   spectra_max_conc=r.df_raw.to_numpy().max())
    r.video = video_link
    # TODO: warn the user if termination is early.
    return r


def scrn_simulate(rxns, time_max=100, lattice=None, display_class=None, video=False, spectra_in_video=True,
                  spectra_average_duration=2, species_tracked=[], manifest_file="", rng_seed=923123122,
                  video_path="", ensemble_size=1):
    video_path = resolve_video(video, video_path)
    if video_path == -1:
        return
    if ensemble_size == 1:
        return scrn_simulate_single_run(rxns, time_max=time_max, lattice=lattice, display_class=display_class,
                                        video=video, spectra_in_video=spectra_in_video,
                                        spectra_average_duration=spectra_average_duration,
                                        species_tracked=species_tracked, manifest_file=manifest_file,
                                        rng_seed=rng_seed, video_path=video_path)
    else:
        ensemble_results = []
        for i in range(ensemble_size):
            run_video_path = f"{video_path}/{i}"
            results = scrn_simulate_single_run(rxns, time_max=time_max, lattice=lattice, display_class=display_class,
                                               video=video, spectra_in_video=spectra_in_video,
                                               spectra_average_duration=spectra_average_duration,
                                               species_tracked=species_tracked, manifest_file=manifest_file,
                                               rng_seed=rng_seed + 2 * i, video_path=run_video_path)
            ensemble_results.append(results)
        return Ensemble(ensemble_results)


def resolve_video(video, video_path):
    if not video:
        return ""
    video_from_argument = True if video_path else False
    if video and not video_from_argument:
        video_path = input(f"Name a directory to store frames and videos: \n{os.getcwd()}/")

        if not video_path:
            video_path = "Surface CRN Videos"
            print(f"Using the default directory {os.getcwd()}/{video_path}")

        # TODO: ask users to press return or enter
        clear_output(wait=False)

    if os.path.isdir(video_path):
        use_path = False
        wrong_decision_word = False

        while not use_path:
            if wrong_decision_word:

                use_path = input(f"Type \"Yes\" to overwrite the directory, or \"No\" if otherwise: ")
                print('\n')
            else:
                use_path = input(f"The directory {os.getcwd()}/{video_path} already exists, would you like to "
                                 f"overwrite the directory? \nType \"Yes\" if you do, or \"No\" if otherwise: ")
                # if same_answer:
                print("\n")

            if use_path.lower() == "yes":
                use_path = True
                wrong_decision_word = False
            elif use_path.lower() == "no":
                if video_from_argument:
                    print("Please choose a different path for the videos.")
                    print("Program exits")
                    return -1

                new_video_path = input(f"Name a directory to store frames and videos: \n{os.getcwd()}/")

                if not new_video_path:
                    new_video_path = "Surface CRN Videos"
                    print(f"Using default directory name {os.getcwd()}/{new_video_path}")

                if os.path.isdir(new_video_path):
                    use_path = False

                # TODO: don't clear if the directory is same as before.
                if new_video_path != video_path:
                    clear_output(wait=False)
                else:
                    same_answer = True
                wrong_decision_word = False
                video_path = new_video_path
            else:
                wrong_decision_word = True
                print(f"\"{use_path}\" is not a valid input.")
                use_path = False
    else:
        os.mkdir(video_path)
    return video_path if video_path else ""


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


def simulate_with_display(manifest_file, lattice, group_selection_seed, rxns=None, spectra_in_video=True,
                          running_average=10, spectra_max_conc=-1):
    if rxns.surface.use_coord_grid:
        display_class = CoordGridDisplay
    elif rxns.surface.structure == "hexagon":
        display_class = HexGridPlusIntersectionDisplay
    else:
        display_class = None
    concs, times = SurfaceCRNQueueSimulator.simulate_surface_crn(manifest_file, group_selection_seed, display_class,
                                                                 init_state=lattice, rxns=rxns,
                                                                 spectra_in_video=spectra_in_video,
                                                                 running_average=running_average,
                                                                 spectra_max_conc=spectra_max_conc)
    return Results.concs_times_df(concs, times)


def simulate_without_display(manifest_file, lattice, species_tracked, rxns, group_selection_seed):
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

    # add_groups(lattice, rxns)
    simulator = QueueSimulator(surface=lattice,
                               transition_rules=opts.transition_rules,
                               seed=opts.rng_seed,
                               group_selection_seed=group_selection_seed,
                               simulation_duration=opts.max_duration,
                               rxns=rxns)

    times = [0]
    concs = dict()
    for species in species_tracked:
        concs[species] = [0]
    for node in lattice:
        if node.state in concs:
            concs[node.state][0] += 1
    ipython_visuals.update_progress(0 / opts.max_duration, "Simulation in progress")

    while not simulator.done():
        # Advance the reaction
        next_rxn = simulator.process_next_reaction()
        if next_rxn is None:
            break
        times.append(next_rxn.time)

        counter = simulator.surface.species_count()
        for species in species_tracked:
            concs[species].append(counter[species])

        # concs[species].append(concs[species][-1])
        # # Very simple mechanism, one reaction at a time.
        # # TODO: these don't account for size 2.
        # for reactant in next_rxn.rule.inputs:
        #     if reactant in concs:
        #         concs[reactant][-1] -= 1
        # for product in next_rxn.rule.outputs:
        #     if product in concs:
        #         concs[product][-1] += 1

        ipython_visuals.update_progress(next_rxn.time/opts.max_duration, "Simulation in progress")

    if times:
        last_time_stamp = times[-1]
    else:
        last_time_stamp = 0
    #     TODO: add rel_tol
    if math.isclose(last_time_stamp, opts.max_duration, abs_tol=0.1):
        ipython_visuals.update_progress(last_time_stamp / opts.max_duration, "Simulation completed", terminating=True)
    else:
        ipython_visuals.update_progress(last_time_stamp / opts.max_duration,
                                        f"Simulation terminated early at {last_time_stamp:.1f} s", terminating=True)
    return times, concs
