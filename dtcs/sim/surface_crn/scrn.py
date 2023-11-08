"""TODO(Andrew)
"""

from dtcs.common import util
from dtcs.sim.surface_crn.surface_crns.simulators.queue_simulator import QueueSimulator
from dtcs.sim.surface_crn.surface_crns.simulators.event_history import EventHistory
from dtcs.sim.surface_crn.surface_crns.base.transition_rule import TransitionRule


def simulate_scrn_single_nodisplay(
        rsys,
        species_names,
        init_surface,
        time_max,
        rng_seed,
        opts,
):
    reactions = [rxn.to_rxns() if hasattr(rxn, 'to_rxns') else rxn for rxn in rsys.reactions]
    transition_rules = []
    for rxn in util.flat(reactions):
        reactants = [react.name for react in rxn.reactants]
        products = [prod.name for prod in rxn.products]
        transition_rules.append(TransitionRule(reactants, products, rxn.get_rate()))

    # Make the simulator
    simulator = QueueSimulator(
        surface=init_surface,
        transition_rules=transition_rules,
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


