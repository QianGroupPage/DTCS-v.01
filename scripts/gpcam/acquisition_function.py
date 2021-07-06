import numpy as np
import matplotlib.pyplot as plt

###############################
#README:
#acquisition functions can be defined here following the templates below
#the acquisition function is a function defined as f=f(x), X --> R
#x is the point where to evaluate the acquisition function
#obj is the gp object that contains functions such as posterior_mean()
#and posterior_covariance(), shannon_information_gain(),...
#the return is either a scalar or a 1d array
#that contains the acquisition function values for a set of points x
#The acquisition funciton defined here will be MAXIMIZED in the algorithm
#to find the next optimal point
###############################

def exploration(x,obj):
    cov = obj.posterior_covariance(x)["v(x)"]
    return cov


def upper_confidence_bounds(x,obj):
    a = 3.0 #####change here, 3.0 for 95 percent confidence interval
    cov = obj.posterior_mean(x)["f(x)"]
    mean = obj.posterior_covariance(x)["v(x)"]
    return mean + a * cov

import argparse
import sys
import os
import json
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import subprocess
import random

from lblcrn.crn_sym import *
from lblcrn.experiments.simulate import simulate
from lblcrn.experiments.xps_io import read_new_data
from lblcrn.experiments.storage import CRNStorage, CRNData

sm = SpeciesManager()

y1 = sm.sp('H2Og', Orbital('1s', 535.0))
x2 = sm.sp('H2O', Orbital('1s', 532.2))
x3 = sm.sp('OH', Orbital('1s', 530.9))
x4 = sm.sp('O', Orbital('1s', 530.0))
x53 = sm.sp('OH.H2O_hb', Orbital('1s', 531.6))
x54 = sm.sp('O.H2O_hb', Orbital('1s', 531.6))
x6 = sm.sp('multiH2O', Orbital('1s', 533.2))
x7 = sm.sp('O2g', Orbital('1s', 535.0))


constants = [3.207654,1.363342,6.220646,0.160755,0.299507,0.167130,1.939313,0.515646,0.733491,0.311754,1.038423, 0.962999,0.002342,426.922895]

def rsys_generator(scaled):
    rsys = RxnSystem(
        Rxn(x4 + y1, x54, scaled[0]),
        Rxn(x3 + y1, x53, scaled[1]),
        Rxn(x54, x3 + x3, scaled[2]),
        Rxn(x3 + x3, x54, scaled[3]),
        Rxn(x53, x2 + x3, scaled[4]),
        Rxn(x54, x2 + x4, scaled[5]),
        Rxn(x2, y1, scaled[6]),
        Rxn(y1, x2, scaled[7]),
        Rxn(x53, y1 + x3, scaled[8]),
        Rxn(x54, x4 + y1, scaled[9]),
        Rxn(x53 + y1, x6, scaled[10]),
        Rxn(x6, x53 + y1, scaled[11]),
        Rxn(x4 + x4, x7, scaled[12]),
        Rxn(x7, x4 + x4, scaled[13]),
        Conc(y1,1),
        Conc(x4,0.25),
        sm
    )
    return rsys

def simulate_and_compare(scaled, exp_env, exp_be):
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

def crn_acquisition():
    data = read_new_data('../../data/1e-1_302k.txt')[0]

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
            xps, rmse = simulate_and_compare(scaled, exp_env, exp_be)
            print(rmse, x[0])
            rmses.append(-rmse)
        return np.array(rmses)

    return acquisition
