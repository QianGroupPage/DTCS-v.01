import numpy as np
import pandas as pd

from lblcrn.twin.core import simulate
from lblcrn.io.xps import read_new_data


def simulate_and_compare(rsys_generator, scaled, exp_env, exp_be):
    rsys = rsys_generator(scaled)
    xps, ts = simulate(rsys, time=500, title="", experimental=pd.Series(data=exp_env, index=exp_be))

    db_env = np.array([])

    raw_db_env_base = xps.df.simulated.envelope.to_numpy()
    raw_db_be_base = xps.df.simulated.index.to_numpy().astype(float)

    raw_db_env, raw_db_be = np.array([]), np.array([])
    for i, be in enumerate(raw_db_be_base):
        if be <= 534: # TODO: Remove this hardcoded constant
            raw_db_env = np.append(raw_db_env, raw_db_env_base[i])
            raw_db_be = np.append(raw_db_be, be)

    input_i = 0

    scale_factor = max(exp_env) / max(raw_db_env)
    raw_db_env *= scale_factor

    # Resample the stored data to match the number of datapoints in the experimental input.
    for j, pt in enumerate(raw_db_be):
        if pt > exp_be[input_i]:
            if j == 0:
                db_env = np.append(db_env, raw_db_env[j])
            else:
                db_env = np.append(db_env, raw_db_env[j-1])
            input_i += 1
        if input_i == len(exp_env):
            break

    if len(db_env) < len(exp_env):
        db_env = np.append(db_env, np.array([db_env[-1] for _ in range(len(exp_env)-len(db_env))]))

    rmse = ((exp_env - db_env)**2).mean() **.5

    return xps, rmse

best_constants_count = 10

class CRNInstrumentation:
    def __init__(self, rsys_generator, experimental_file_path):
        data = read_new_data(experimental_file_path)[0]

        intensities, bes = np.array([]), np.array([])
        for i, be in enumerate(data.binding_energies):
            if be <= 534:
                intensities = np.append(intensities, data.intensities[i])
                bes = np.append(bes, be)

        series = pd.Series(data=intensities, index=bes)

        self.exp_env = series.to_numpy()
        self.exp_be = series.index.to_numpy().astype(float)
        self.rsys_generator = rsys_generator

        self.rmse_evolution = []
        # Contains tuples of (constants, rmse).
        # This list will be ordered in ascending order of rmse.
        self.best_constants = []

    def func(self):
        def instrumentation(data):
            calculated_rmses = []
            for instance in data:
                scaled = instance["position"]
                xps, rmse = simulate_and_compare(self.rsys_generator, scaled, self.exp_env, self.exp_be)

                if len(self.best_constants) < best_constants_count:
                    self.best_constants.append((list(scaled), rmse))
                    self.best_constants.sort(key=lambda x: x[1])
                elif max(self.best_constants, key=lambda x: x[1])[1] > rmse:
                    self.best_constants[len(self.best_constants)-1] = (list(scaled), rmse)
                    self.best_constants.sort(key=lambda x: x[1])

                calculated_rmses.append(rmse)
                instance["value"] = rmse
                print(rmse)

            self.rmse_evolution.append(np.min(calculated_rmses))

            return data
        return instrumentation
