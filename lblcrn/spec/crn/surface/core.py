"""TODO"""

from typing import Optional

import random

from lblcrn.spec.crn.crn_abc import CRNSpecABC
from lblcrn.spec.crn.surface.rxn_system import SurfaceRxnSystem
from lblcrn.spec.crn.surface.species import SurfaceSpeciesManager


class SurfaceCRNSpec(CRNSpecABC):
    """TODO
    """

    _rxn_sys_cls = SurfaceRxnSystem

    def __init__(self,
                 *components,
                 rsys: SurfaceRxnSystem = None,
                 species: SurfaceSpeciesManager = None,
                 time: int = 10,
                 runs: int = 1,
                 #rng_seed: Optional[int] = None,
                 **kwargs):
        super().__init__(*components,
                         rsys=rsys,
                         species=species,
                         sim_type='surface')
        self.time = time
        self.runs = runs

        # rng_seed = rng_seed or random.randrange(1024)
        #
        # manifest = generate_manifest_stream(
        #     rsys=rsys,
        #     max_duration=self.time,
        #     random_seed_scrn=rng_seed,
        #     random_seed_surface=rng_seed,
        #     video_path='video path'
        # )
        #
        # manifest_opts = read_manifest(
        #     filename=manifest,
        # )
        #
        # opts = SurfaceCRNOptionParser(
        #     options=manifest_opts
        # )
        #
        # self.__dict__.update(opts.__dict__)

        # Ye's
        # self.movie_title = 'SCRN Simulation'
        # self.speedup_factor = 0.5
        # self.debug = False
        # self.fps = 1
        # self.display_text = True
        # self.pixels_per_node = 80
        # self.wrap_grid = False
        # self.rng_seed = random.randrange(1024)
        # self.max_duration = time
        # self.capture_directory = 'video path'
        #
        # self.capture_rate
        # self.COLORMAP
        # self.simulation_type
        # if self.simulation_type
        #     self.update_rule
        # elif self.simulation_type
        #     self.transition_rules
        # self.surface_geometry
        # self.grid_type
        # if self.grid_type
        #     self.emulation_colormap
        #     self.horizontal_buffer
        #     self.vertical_buffer
        #     self.cell_height
        #     self.cell_width
        #     self.representative_cell_x
        #     self.representative_cell_y
        # self.init_state
