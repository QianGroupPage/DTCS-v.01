import time

from gpcam.autonomous_experimenter import AutonomousExperimenterGP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lblcrn.bulk_crn.gpcam.acquisition import CRNAcquisition
from lblcrn.bulk_crn.gpcam.kernel import kernel_l2_single_task

def evaluate(rsys_generator, constants, experimental_file_path):
    acquisition = CRNAcquisition(rsys_generator, experimental_file_path)

    exp = AutonomousExperimenterGP(
        parameter_bounds=np.array([[0.1, 10.0] for _ in range(len(constants))]),
        instrument_func=synthetic_instrumentation_function,
        hyperparameters=[1.0 for _ in range(1+len(constants))],
        hyperparameter_bounds=[[1.0,100.0] for _ in range(1+len(constants))],
        init_dataset_size=len(constants),
        acq_func=acquisition.func(),
        kernel_func=kernel_l2_single_task, # TODO(rithvikp): Add configurability.
    )

    exp.go(N=20)

    return acquisition

# Credit to the GPCam package for the following functions.
def synthetic_instrumentation_function(data):
    for idx_data in range(len(data)):
        if data[idx_data]["measured"] == True: continue
        x1 = data[idx_data]["position"][0]
        x2 = data[idx_data]["position"][1]
        data[idx_data]["value"] = himmel_blau([x1, x2])
        data[idx_data]["cost"] = {"origin": np.random.uniform(low=0.0, high=1.0, size = 2),
                                     "point": np.random.uniform(low=0.0, high=1.0, size = 2),
                                     "cost": np.random.uniform(low=0.0, high=1.0)}
        data[idx_data]["variance"] = 0.01
        data[idx_data]["measured"] = True
        data[idx_data]["time stamp"] = time.time()
    return data

def himmel_blau(x):
    return (x[0] ** 2 + x[1] - 11.0) ** 2 + (x[0] + x[1] ** 2 - 7.0) ** 2
