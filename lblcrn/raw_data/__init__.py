"""A module for reading and processing raw data from XPS experiments.
"""

from lblcrn.raw_data.raw_experiment import RawExperiment
from lblcrn.raw_data.raw_region import RawRegion
from lblcrn.raw_data.raw_measurement import RawMeasurement
from lblcrn.raw_data.decomposition import decompose
from lblcrn.raw_data.baseline import shirley_background
