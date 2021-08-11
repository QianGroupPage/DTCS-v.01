import time

from gpcam.autonomous_experimenter import AutonomousExperimenterGP
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lblcrn.bulk_crn.gpcam.instrumentation import CRNInstrumentation
from lblcrn.bulk_crn.gpcam.kernel import kernel_l2_single_task

def evaluate(rsys_generator, constants, experimental_file_path):
    instrumentation = CRNInstrumentation(rsys_generator, experimental_file_path)

    exp = AutonomousExperimenterGP(
        parameter_bounds=np.array([[0.1, 10.0] for _ in range(len(constants))]),
        instrument_func=instrumentation.func(),
        hyperparameters=np.ones((1+len(constants),)),
        hyperparameter_bounds=[[0.1,100.0] for _ in range(1+len(constants))],
        #acq_func = "covariance",
        init_dataset_size=len(constants),
        kernel_func=kernel_l2_single_task, # TODO(rithvikp): Add configurability.
    )

    exp.train()
    exp.go(N=25)

    return instrumentation
