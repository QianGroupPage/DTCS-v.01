"""TODO"""

#GENERAL IMPORTS
# from dtcs import *
# from gpcam import *

#INSTRUMENTATION FUNCTION IMPORTS
# import argparse
# from datetime import datetime, timedelta
# import json
# import os
# import random
# import subprocess
# import sys
import scipy

# from dtcs.io.xps import read_new_data
# from dtcs.io.xps import read_new_data
# from dtcs.io.storage import CRNStorage, CRNData

#EVALUATION FUNCTION
from gpcam.autonomous_experimenter import AutonomousExperimenterGP

# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#MAIN TEST IMPORTS
# from dtcs.spec.crn.rxn_abc import *
# from dtcs.spec.crn.bulk.rxn_system import *
# from dtcs.spec.species import *
# from dtcs.spec.xps import *
# from dtcs.sim.bulk_crn.bulk_crn import *
# from dtcs.spec.crn.bulk.reaction import BulkRxn as Rxn
# from dtcs.spec.crn.bulk.core import BulkCRNSpec


class CRNInstrumentation:
    def __init__(self, crn, experimental, ignore):
        # data = read_new_data(experimental_file_path)[0]
        #
        # intensities, bes = np.array([]), np.array([])
        # for i, be in enumerate(data.binding_energies):
        #     if be <= 540:
        #         intensities = np.append(intensities, data.intensities[i])
        #         bes = np.append(bes, be)
        #
        # series = pd.Series(data=intensities, index=bes)

        self.exp_env = experimental.to_numpy()
        self.exp_be = experimental.index.to_numpy().astype(float)
        self.crn = crn
        self.rsys_generator = crn.subs_rates
        self.ignore = ignore

        self.rmse_evolution = []
        # Contains tuples of (constants, rmse).
        # This list will be ordered in ascending order of rmse.
        self.best_constants = []
        # good constants
        self.constants = []
        # bulk list of simulated xps's I think
        self.xpss = []

        # user information
        self.bestr2 = -1
        self.bestxps = None
        self.location = 0

    def func(self):
        def instrumentation(data):
            calculated_rmses = []
            c=0
            if(len(data) < 15):
                for instance in data:
                    vals = instance['position']

                    # scaled = [vals[0], vals[1], vals[2], 1/vals[2], vals[3], vals[4], vals[5], 1/vals[5], 1/vals[1], 1/vals[0], vals[6], 1/vals[6], vals[7], 1/vals[7]]
                    scaled = vals

                    xps, r_squared, rmse = simulate_and_compare(
                        crn=self.crn,
                        scaled=scaled,
                        exp_env=self.exp_env,
                        exp_be=self.exp_be,
                        ignore=self.ignore
                    )
                    instance['value'] = r_squared

                    if(r_squared > self.bestr2):
                        self.bestr2 = r_squared
                        self.bestxps = xps
                        self.constants = scaled
                        self.location=c

                xps, r_squared, rmse = simulate_and_compare(
                    crn=self.crn,
                    scaled=scaled,
                    exp_env=self.exp_env,
                    exp_be=self.exp_be,
                    ignore=self.ignore
                )

                print('RMSE: ' + str(rmse) + ', ' + 'r^2: ' + str(r_squared))
                # self.rmse_evolution.append(min(calculated_rmses))
            else:
                instance = data[len(data)-1]
                vals = instance['position']
                #           0         1        2        3          4        5        6       7          8            9       10          11         12        13
                # scaled = [vals[0], vals[1], vals[2], 1/vals[2], vals[3], vals[4], vals[5], 1/vals[5], 1/vals[1], 1/vals[0], vals[6], 1/vals[6], vals[7], 1/vals[7]]
                scaled = vals

                xps, r_squared, rmse = simulate_and_compare(
                    crn=self.crn,
                    scaled=scaled,
                    exp_env=self.exp_env,
                    exp_be=self.exp_be,
                    ignore=self.ignore
                )
                instance['value'] = r_squared

                if(r_squared > self.bestr2):
                    self.bestr2 = r_squared
                    self.bestxps = xps
                    self.constants = scaled
                    self.location = c

                print('RMSE: ' + str(rmse) + ', ' + 'r^2: ' + str(r_squared))
                # self.rmse_evolution.append(min(calculated_rmses))


            return data
        return instrumentation


def kernel_l2_single_task(x1, x2, hyperparameters, obj):
    ################################################################
    ###standard anisotropic kernel in an input space with l2########
    ###########################for single task######################
    '''
    x1: 2d numpy array of points
    x2: 2d numpy array of points
    obj: object containing kernel definition

    Return:
    -------
    Kernel Matrix
    '''
    hps = hyperparameters
    distance_matrix = np.zeros((len(x1), len(x2)))
    # you can always get some help:
    # help(obj.squared_exponential_kernel)
    for i in range(len(x1[0]) - 1):
        distance_matrix += abs(np.subtract.outer(x1[:, i], x2[:, i]) / hps[1 + i]) ** 2
    distance_matrix = np.sqrt(distance_matrix)
    # return   hps[0] * obj.squared_exponential_kernel(distance_matrix,1)
    # return   hps[0] * obj.exponential_kernel(distance_matrix,1) + obj.periodic_kernel(abs(x1[:,-1]-x2[:,-1]), hps[-1],hps[-2])
    return hps[0] * obj.matern_kernel_diff1(distance_matrix, 1)
    # return   hps[0] * obj.matern_kernel_diff2(distance_matrix,1)


def simulate_and_compare(crn, scaled, exp_env, exp_be, ignore):
    crn = crn.subs_rates(scaled)
    xps = crn.simulate_xps(title='',
                           experimental=pd.Series(data=exp_env, index=exp_be),
                           ignore=ignore)

    sim_env = np.array([])
    reconv_env = np.array([])

    raw_reconv_env_base = xps.df.deconvoluted.envelope.to_numpy()
    raw_reconv_be_base = xps.df.deconvoluted.index.to_numpy().astype(float)

    raw_sim_env_base = xps.df.simulated.envelope.to_numpy()
    raw_sim_be_base = xps.df.simulated.index.to_numpy().astype(float)

    raw_sim_env, raw_sim_be = np.array([]), np.array([])
    raw_reconv_env, raw_reconv_be = np.array([]), np.array([])

    for i, be in enumerate(raw_sim_be_base):
        if be <= 540:
            raw_sim_env = np.append(raw_sim_env, raw_sim_env_base[i])
            raw_reconv_env = np.append(raw_reconv_env, raw_reconv_env_base[i])
            raw_sim_be = np.append(raw_sim_be, be)
            raw_reconv_be = np.append(raw_reconv_be, be)

    input_i = 0

    # Resample the stored data to match the number of datapoints in the experimental input.
    for j, pt in enumerate(raw_sim_be):
        if pt > exp_be[input_i]:
            if j == 0:
                sim_env = np.append(sim_env, raw_sim_env[j])
                reconv_env = np.append(reconv_env, raw_reconv_env[j])
            else:
                sim_env = np.append(sim_env, raw_sim_env[j - 1])
                reconv_env = np.append(reconv_env, raw_reconv_env[j - 1])
            input_i += 1
        if input_i == len(exp_env):
            break

    if len(sim_env) < len(exp_env):
        sim_env = np.append(sim_env, np.array([sim_env[-1] for _ in range(len(exp_env) - len(sim_env))]))
        reconv_env = np.append(reconv_env, np.array([reconv_env[-1] for _ in range(len(exp_env) - len(reconv_env))]))

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(reconv_env, sim_env)

    rmse = ((exp_env - sim_env) ** 2).mean() ** .5
    r_squared = r_value ** 2

    return xps, r_squared, rmse


def evaluate(crn, experimental, iterations, ignore):
    num_rates = len(crn.get_rates())
    instrumentation = CRNInstrumentation(
        crn=crn,
        experimental=experimental,
        ignore=ignore,
    )

    hpbounds = [[0.001, 1000]] + [[0.1, 5000] for _ in range(num_rates)]
    exp = AutonomousExperimenterGP(
        parameter_bounds=np.array([[0.01, 100] for _ in range(num_rates)]),
        instrument_func=instrumentation.func(),
        hyperparameters=np.ones(num_rates + 1),
        hyperparameter_bounds=hpbounds,
        acq_func='ucb',
        init_dataset_size=num_rates,
        kernel_func=kernel_l2_single_task
    )

    exp.train()
    exp.go(
        N=iterations,
        acq_func_opt_setting=lambda number: 'global' if number % 2 == 0 else 'global',
   )

    return instrumentation
