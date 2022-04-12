from __future__ import annotations

import copy
import warnings

from dtcs.sim import surface_crn
from dtcs.spec.crn.surface.api_adapt import generate_surface, generate_manifest_stream
from dtcs.sim.surface_crn.surface_crns.options.option_processor import SurfaceCRNOptionParser
from dtcs.sim.surface_crn.surface_crns.readers.manifest_readers import read_manifest
from dtcs.sim.surface_crn.surface_crns.simulators.event_history import EventHistory
from dtcs.sim.surface_crn.surface_crns.simulators.queue_simulator import QueueSimulator
from dtcs.twin.crn import SurfaceCRNTimeSeries


def depr_simulate_surface_crn(scrn, **kwargs):
    ens = surface_crn.scrn.scrn_simulate(scrn.rsys,
                                         time_max=scrn.time,
                                         ensemble_size=scrn.runs,

                                         video=False,
                                         spectra_in_video=False,
                                         video_path='output',
                                         trajectory_path='output',)
    warnings.warn('TODO', DeprecationWarning)
    results = None
    if scrn.runs == 1:
        results = [ens]
    else:
        results = ens.results

    scts = SurfaceCRNTimeSeries.from_runs([result.df_raw for result in ens.results], scrn)
    return scts


def _simulate_surface_crn(
        scrn: SurfaceCRNSpec,
        *,
        rsys: SurfaceRxnSystem = None,
        species: SurfaceSpeciesManager = None,

        surface: Surface = None,

        time: float = None,
        num_runs: int = None,
        initial_state_seed: int = None, # TODO(Andrew) should there be separate seeds?
        # the initial state should be different across different runs, right?
        simulation_rng_seed: int = None,
):

    # --- Extract input from Spec -----------------------------------------
    rsys = rsys or scrn.rsys
    species = species or scrn.species
    surface = surface or scrn.surface
    # TODO(Andrew) Allow for initial surface state (see if it's in the Surface class)

    # Collect the names of all the species
    # TODO(Andrew): the SurfaceSpeciesManager should be know about the surface
    species_names = copy.copy(scrn.species.names)
    species_names.append(surface.name)
    species_names.extend([site.name for site in surface.sites])

    time = time or scrn.time
    num_runs = num_runs or scrn.runs # TODO(Andrew) change name

    initial_state_seed = initial_state_seed or scrn.initial_state_seed
    simulation_rng_seed = simulation_rng_seed or scrn.simulation_rng_seed

    # --- Run num_runs simulations ------------------------------------
    runs_init_surfaces = []
    runs_times = []
    runs_species_counts = []
    event_histories = []

    for run_number in range(num_runs):
        run_rng_seed = simulation_rng_seed + run_number  # TODO(Andrew) this might be bad if someone runs several scrns

        # Generate an initial surface randomly if one isn't supplied
        init_surface_state = generate_surface(rsys)  # TODO(Andrew) this function only needs rsys.surface
        # TODO(Andrew) it also might not even be used (???)
        # TODO(Andrew) use initial_state_seed here; this should be random (?)
        runs_init_surfaces.append(copy.deepcopy(init_surface_state))  # TODO(Andrew) best way to track this?

        times, species_counts, event_history = _simulate_scrn_single_nodisplay(
            rsys=rsys,
            species_names=species_names,
            init_surface=init_surface_state,
            time_max=time,
            rng_seed=run_rng_seed
        )

        runs_times.append(times)
        runs_species_counts.append(species_counts)
        event_histories.append(event_history)

    scts = SurfaceCRNTimeSeries.from_times_and_counts(runs_times, runs_species_counts, scrn)
    # TODO(Andrew) make this part of the actual SurfaceCRNTimeSeries.__init__:
    scts._event_histories = event_histories
    scts._init_surfaces = runs_init_surfaces
    return scts


def _simulate_scrn_single_nodisplay(
        rsys,
        species_names,
        init_surface,
        time_max,
        rng_seed,
):
    # Convert input into SurfaceCRNOptionParser format for simulator (sclamons)
    # TODO(Andrew) this is temporary, I should dissect it
    #  as of current, if something is available outside of the manifest, I try to
    #  get it from there. This should ultimately be a method of SCRNSpec
    opts = _get_opts_via_manifest(
        rsys=rsys,
        time=time_max,
        rng_seed=rng_seed
    )

    # Make the simulator
    simulator = QueueSimulator(
        surface=init_surface,
        transition_rules=opts.transition_rules,
        seed=rng_seed,
        group_selection_seed=rng_seed * 2,  # TODO(Andrew) how many different seed options should I give?
        simulation_duration=time_max,
        rxns=rsys,
    )

    # Prepare to store data while simulating
    def count_species():
        counter = simulator.surface.species_count()
        return {specie: counter[specie] for specie in species_names}

    times = [0]
    species_counts = [count_species()]
    event_history = EventHistory()

    # --- Simulate -----------------
    while not simulator.done():
        # Get the next reaction
        next_reaction = simulator.process_next_reaction()
        if next_reaction is None:
            break

        # Record: the time each reaction happens
        times.append(next_reaction.time)

        # Record: the count of each species at each time
        species_counts.append(count_species())

        # Record: each reaction as it happens
        event_history.add_event(next_reaction)
        event_history.increment_event(1)

    return times, species_counts, event_history.history


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

    return opts  #.__dict__