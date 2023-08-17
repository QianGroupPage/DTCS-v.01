"""TODO"""

from collections import defaultdict
import itertools
import numpy as np

# A utilty function to remove unnecessary information
def prune_row_dict(dic):
    first_value = dic[next(iter(dic))]

    # Base case: actual data can have just one column
    if not isinstance(first_value, dict):
        return dic

    # If the dictionary is one-element, drop that row.
    if len(dic) == 1:
        return prune_row_dict(first_value)

    # Otherwise, recurse
    return {key: prune_row_dict(value) for key, value in dic.items()}


def system_generator_conc(
        crn,
        sm,
        gibbs_to_rates,
        temperatures=(298,),  # K
        pressures=(0.1,),  # Torr
        times=(-1,),
        conds=None,  # This is of the format (temp, pressure, time) if you don't want the Cartesian product
        prune=True,
):
    # Sanitize times
    times = tuple(time if time > 0 else crn.time for time in times)
    # If conditions are supplied exactly, don't prune.
    if conds:
        prune = False
    # Create combintions
    conds = conds or tuple(itertools.product(temperatures, pressures, times))

    def rescale_concs(raw):
        peak_dict = defaultdict(float)
        for specie, conc in raw.items():
            if specie.name == 'H2Og':
                continue
            for orbital in sm[specie].orbitals:
                peak_dict[orbital.binding_energy] += conc * orbital.splitting
        return raw / max(peak_dict.values())

    def _sim_at_conds(crn, gibbs, temp, pressure, times: tuple):
        rates = gibbs_to_rates(
            gibbs=gibbs,
            temp=temp,
            pressure=pressure,
        )

        cts_g = crn.subs_rates(rates=rates).simulate(
            time=max(times + (crn.time,))
        )
        # TODO: This rescaling is not well-integrated
        scaled_concs = dict()
        for time in times:
            scaled_concs[time] = rescale_concs(cts_g.at(time))
        return scaled_concs

    def system(gibbs):
        # Simulate the CRN
        out_dict = defaultdict(lambda: defaultdict(dict))
        for (temp, pressure), group in itertools.groupby(conds, key=lambda cond: (cond[0], cond[1])):
            times = tuple(time for _, _, time in group)
            concs_raw = _sim_at_conds(
                crn=crn,
                gibbs=gibbs,
                temp=temp,
                pressure=pressure,
                times=times,
            )

            # Sanitize the output data
            for time in times:
                concs = {symbol.name: conc for symbol, conc in concs_raw[time].items()}
                # concs['O(H)-H2O'] = concs['OH-H2O'] + concs['O-H2O']
                out_dict[f'{temp} K'][f'{pressure} Torr'][f'{time} s'] = concs

        if prune:
            return prune_row_dict(out_dict)
        else:
            return out_dict

    return system


def system_generator(
        crn,
        # gibbs_to_rates,
        sample_at_ev,
        temperatures=(298,),  # K
        pressures=(0.1,),  # Torr
        times=(-1,),
        conds=None,  # This is of the format (temp, pressure, time) if you don't want the Cartesian product
        prune=True,
):
    sample_at_ev = np.array(sample_at_ev)
    sm = crn.sm

    # Sanitize times
    times = tuple(time if time > 0 else crn.time for time in times)
    # If conditions are supplied exactly, don't prune.
    if conds:
        prune = False
    # Create combintions
    conds = conds or tuple(itertools.product(temperatures, pressures, times))

    def rescale_samples(raw):
        return raw / np.max(raw)

    def _sim_at_conds(crn, gibbs, temp, pressure, times: tuple):
        cts_g = crn.subs_gibbs(gibbs=gibbs).simulate(
            time=max(times + (crn.time,)),
            temperature=temp,
            pressure=pressure,
        )

        # TODO: This rescaling is not well-integrated
        scaled_samples = dict()
        for time in times:
            xo = cts_g.xps_with(t=time, autoresample=False).resample(x_range=sample_at_ev)

            scaled_samples[time] = rescale_samples(xo.sim_envelope)
        return scaled_samples

    def system(gibbs):
        # Simulate the CRN
        out_dict = defaultdict(lambda: defaultdict(dict))
        for (temp, pressure), group in itertools.groupby(conds, key=lambda cond: (cond[0], cond[1])):
            times = tuple(time for _, _, time in group)
            samples_raw = _sim_at_conds(
                crn=crn,
                gibbs=gibbs,
                temp=temp,
                pressure=pressure,
                times=times,
            )

            # Sanitize the output data
            for time in times:
                out_dict[f'{temp} K'][f'{pressure} Torr'][f'{time} s'] = dict(samples_raw[time])

        if prune:
            return prune_row_dict(out_dict)
        else:
            return out_dict

    return system

def xps_generator(
        crn,
        temperatures=(298,),  # K
        pressures=(0.1,),  # Torr
        times=(-1,),
        conds=None,  # This is of the format (temp, pressure, time) if you don't want the Cartesian product
        prune=True,
):
    """TODO: in a perfect world, this would return a collection of XPSSpecs."""
    # Sanitize times
    times = tuple(time if time > 0 else crn.time for time in times)
    # If conditions are supplied exactly, don't prune.
    if conds:
        prune = False
    # Create combintions
    conds = conds or tuple(itertools.product(temperatures, pressures, times))

    def _xps_at_conds(crn, gibbs, temp, pressure, times: tuple):
        cts_g = crn.subs_gibbs(gibbs=gibbs).simulate(
            time=max(times + (crn.time,)),
            temperature=temp,
            pressure=pressure,
        )

        xps_by_time = dict()
        for time in times:
            xps_by_time[time] = cts_g.xps_with(t=time, autoresample=False)
        return xps_by_time

    def system(gibbs):
        # Simulate the CRN
        out_dict = defaultdict(lambda: defaultdict(dict))
        for (temp, pressure), group in itertools.groupby(conds, key=lambda cond: (cond[0], cond[1])):
            times = tuple(time for _, _, time in group)
            xps_by_time = _xps_at_conds(
                crn=crn,
                gibbs=gibbs,
                temp=temp,
                pressure=pressure,
                times=times,
            )

            # Sanitize the output data
            for time in times:
                out_dict[f'{temp} K'][f'{pressure} Torr'][f'{time} s'] = xps_by_time[time]

        if prune:
            return prune_row_dict(out_dict)
        else:
            return out_dict

    return system