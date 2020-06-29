"""Core functionality for simulating chemistry experiments.

Exports:
    solve_rsys_ode(): Gives a time-series solution for a RxnSystem by solving
        its ODE-system.

Usage:
    Mostly internal.
"""

import collections

import numpy as np
from scipy import integrate


def solve_rsys_ode(rsys, time_max: float = 1, **options):
    """Simulate the given reaction system over time.

    Private in case we want to add multiple possible solving methods.

    Args:
        rsys: ReactionsSystem, the reaction system to simulate
        time_max: The time until which to simulate.
        **options: Forwarded to scipy.integrate.solve_ivp

    Returns:
        A Solution object describing the solution.
    """

    ode_func = rsys.get_ode_functions()
    num_species = len(rsys._symbols)

    # schedule, a dictionary {time : [amount to add for species no. index]}
    schedule = collections.defaultdict(lambda: [0] * num_species)
    for index in range(num_species):
        for time, amount in rsys.scheduler[index].items():
            schedule[time][index] += amount

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
        for index, amount in enumerate(schedule[current_time]):
            current_concs[index] += amount

        partial_sol = integrate.solve_ivp(ode_func, (current_time, next_time),
                                          current_concs, **options)

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
    for index, func in rsys.get_conc_functions().items():
        for tindex in range(sol_t.size):
            sol_y[index][tindex] = func(sol_t[tindex], sol_y[:, tindex])

    return sol_t, sol_y
