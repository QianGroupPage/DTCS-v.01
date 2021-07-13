from gpcam.autonomous_experimenter import AutonomousExperimenterGP
import matplotlib.pyplot as plt
import numpy as np

from acquisition_function import crn_acquisition
from instrument_function import synthetic_function

num_equations = 14

def run():
    exp = AutonomousExperimenterGP(
        parameter_bounds=np.array([[0.1, 10] for _ in range(num_equations)]),
        instrument_func=synthetic_function,
        hyperparameters=[1.0, 1.0, 1.0],
        hyperparameter_bounds=[[1.0,100.0],[0.10,100.0],[0.10,100.0]],
        init_dataset_size=num_equations,
        acq_func=crn_acquisition(),
    )

    optimized = exp.go(N=20)

run()
