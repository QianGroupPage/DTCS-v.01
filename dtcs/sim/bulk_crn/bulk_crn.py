"""Core functionality for simulating chemistry experiments.

Exports:
    solve_rsys_ode(): Gives a time-series solution for a RxnSystem by solving
        its ODE-system.

Usage:
    Mostly internal.
"""

from typing import Callable, Dict, List, TypeAlias
import collections

import numpy as np
from scipy import integrate


Time: TypeAlias = float


def solve_rsys_ode(
        evolution_func: Callable,
        concentration_func: Callable,
        schedule: Dict[Time, Dict[str, float]],
        species: List[str],
        time_max: float,
        end_when_settled: bool,
        method: str='LSODA',
        **options):
    """Simulate the given reaction system over time.

    Private in case we want to add multiple possible solving methods.

    Args:
        rsys: ReactionsSystem, the reaction system to simulate
        time_max: The time until which to simulate.
        **options: Forwarded to scipy.integrate.solve_ivp

    Returns:
        A Solution object describing the solution.
    """

    num_species = len(species)
    # schedule, a dictionary {time : [amount to add for species no. index]}
    # schedule = collections.defaultdict(lambda: [0] * num_species)
    # for index in range(num_species):
    #     for time, amount in rsys.scheduler[index].unpack():
    #         schedule[time][index] += amount

    # This is an ordered list of all the times at which we add/remove stuff.
    time_breaks = sorted(schedule.keys())

    # Do the simulation in broken pieces
    current_concs = [0] * num_species
    current_time = time_breaks.pop(0)  # Guaranteed to be 0

    while current_time < time_max:
        # Get next_time, the end of the interval we're simulating this step.
        if len(time_breaks) > 0:
            next_time = time_breaks.pop(0)
        else:
            next_time = time_max

        # Add the respective amounts for this timestep.
        for index, specie in enumerate(species):
            current_concs[index] += schedule[current_time][specie]

        events = None
        if end_when_settled:
            events = settle_event_creator()

        partial_sol = integrate.solve_ivp(evolution_func,  (current_time, next_time),
                                          current_concs, method=method, events=events, **options)

        # Add the partial solution to the whole solution.
        if current_time == 0:
            sol_t = partial_sol.t
            sol_y = partial_sol.y
        else:
            sol_t = np.append(sol_t, partial_sol.t)
            sol_y = np.append(sol_y, partial_sol.y, axis=1)

        # Loop; set current_concs to the new ones and curren_time to the next one
        for index in range(num_species):
            current_concs[index] = sol_y[index][len(sol_y[index]) - 1]
        current_time = next_time

    # Set y for concentration equations, which the ODE solver does't calculate.
    for index in range(sol_y.shape[1]):
        sol_y[:, index] = concentration_func(sol_t[index], sol_y[:, index])

    return sol_t, sol_y


SETTLE_EVENT_THRESHOLD = 5e-7


def settle_event_creator():
    last_y = None
    
    def settle(t, y):
        nonlocal last_y
        if last_y is None:
            last_y = np.array(y)
        diff = abs(np.sum(np.array(y) - last_y))
        if diff < SETTLE_EVENT_THRESHOLD:
            return 0
        last_y = np.array(y)
        print(diff)
        return diff
    settle.terminal = True
    return settle


def solve_rsys_ode_till_eq(rsys, **options):
    return solve_rsys_ode(rsys, time_max=100, **options)
