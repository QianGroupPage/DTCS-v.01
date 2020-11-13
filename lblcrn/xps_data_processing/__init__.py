"""A module for reading and processing raw data from XPS experiments.
"""

from lblcrn.xps_data_processing.baseline import shirley_background
from lblcrn.xps_data_processing.fitting_suggestions import \
    suggest_fitting_params
from lblcrn.xps_data_processing.peak_fit import decompose
from lblcrn.xps_data_processing.raw_experiment import RawExperiment
from lblcrn.xps_data_processing.raw_measurement import RawMeasurement
from lblcrn.xps_data_processing.raw_region import RawRegion
