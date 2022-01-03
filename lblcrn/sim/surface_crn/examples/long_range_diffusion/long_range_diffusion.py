from surface_crns import SurfaceCRNQueueSimulator
from surface_crns.models.grids import *
from surface_crns.base.node import *
from surface_crns.views.grid_display import *
from surface_crns.simulators.queue_simulator import *
import surface_crns.readers as readers
import pygame, math, random, sys
import matplotlib.pyplot as plt

import lblcrn.sim.surface_crn.core


def main():
    manifest_file = "long_range_diffusion.txt"
    manifest_options = \
                readers.manifest_readers.read_manifest(manifest_file)
    init_state = manifest_options["init_state"]

    lattice = SquareGridWithCornerLeak(init_state.shape[0], init_state.shape[1],
                                       0.25)
    lattice.set_global_state(init_state)
    lblcrn.sim.surface_crn.core._simulate_surface_crn(manifest_file,
                                                      init_state = lattice)

if __name__ == "__main__":
    main()