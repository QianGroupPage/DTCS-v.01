import argparse
from datetime import datetime, timedelta
import json
import os
import random
import subprocess
import sys

from gpcam.autonomous_experimenter import AutonomousExperimenterGP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lblcrn.experiments.simulate import simulate
from lblcrn.experiments.xps_io import read_new_data
from lblcrn.experiments.storage import CRNStorage, CRNData

def simulate_and_compare(rsys_generator, scaled, exp_env, exp_be):
    rsys = rsys_generator(scaled)
    xps, ts = simulate(rsys, time=500, title="")

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

def crn_acquisition(rsys_generator, experimental_file_path):
    data = read_new_data(experimental_file_path)[0]

    intensities, bes = np.array([]), np.array([])
    for i, be in enumerate(data.binding_energies):
        if be <= 534:
            intensities = np.append(intensities, data.intensities[i])
            bes = np.append(bes, be)

    series = pd.Series(data=intensities, index=bes)
    exp_env = series.to_numpy()
    exp_be = series.index.to_numpy().astype(float)

    def acquisition(x, obj):
        rmses = []
        for scaled in x:
            xps, rmse = simulate_and_compare(rsys_generator, scaled, exp_env, exp_be)
            print(rmse, x[0])
            rmses.append(-rmse)
        return np.array(rmses)

    return acquisition
