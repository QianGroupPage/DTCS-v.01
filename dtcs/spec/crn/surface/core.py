"""TODO"""

from __future__ import annotations

import collections
import copy

from dtcs.common.display import color_map
from dtcs.spec.crn.surface.api_adapt import generate_manifest_stream, \
    _get_opts_via_manifest
from dtcs.spec.crn.surface.conditions import Coverage
from dtcs.spec.crn.bulk import CRNSpec, Conc
from dtcs.spec.crn.crn_abc import CRNSpecABC
from dtcs.spec.crn.surface.rxn_system import SurfaceRxnSystem
from dtcs.sim.surface_crn.surface_crns.readers.manifest_readers import read_manifest
from dtcs.sim.surface_crn.surface_crns.options.option_processor import SurfaceCRNOptionParser
from dtcs.sim.surface_crn.surface_crns.models.grids import SquareGrid
from dtcs.sim.surface_crn.hex_grid_with_intersect import HexGridPlusIntersections
from dtcs.sim.surface_crn.surface_crns.models.coord_grid import CoordGrid
from dtcs.sim.surface_crn.scrn import simulate_scrn_single_nodisplay

from dtcs.twin.crn import SurfaceCRNTimeSeries


class SurfaceCRNSpec(CRNSpecABC):
    """TODO
    """

    _rxn_sys_cls = SurfaceRxnSystem

    def __init__(self,
                 *components,
                 rsys: SurfaceRxnSystem = None,
                 species: SpeciesManager = None,
                 surface: Surface = None,
                 size=None,
                 time: int = 10,
                 runs: int = 1,
                 #rng_seed: Optional[int] = None,
                 **kwargs):
        super().__init__(*components,
                         rsys=rsys,
                         species=species,)
        self.time = time
        self.runs = runs
        self.size = size
        self.surface = surface

        # TODO(Andrew) drop the hardcoded numbers
        #  and find a reasonable way to input seeds.
        self.initial_state_seed = 10
        self.simulation_rng_seed = 20

        # TODO(Andrew): Bodges
        self.rsys.surface = surface
        self.surface.size = size or (10, 10)

        # TODO(Andrew): Use _get_opts_via_manifest?
        manifest = generate_manifest_stream(
            rsys=self.rsys,
            max_duration=self.time,
            random_seed_scrn=self.simulation_rng_seed,
            random_seed_surface=self.initial_state_seed,
            video_path='video path'
        )

        manifest_opts = read_manifest(
            filename=manifest,
        )

        opts = SurfaceCRNOptionParser(
            options=manifest_opts
        )

        self.__dict__.update(opts.__dict__)

        # TODO(Andrew) Change the name to not be all caps, jeez
        self.COLORMAP = {species: color_map.rgb256(species) for species in
                         self.sm.names}
        self.COLORMAP.update({site: color_map.rgb256(site) for site
                              in self.surface.sites})

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

    # dtcs.twin.crn
    def simulate(
            self: SurfaceCRNSpec,
            *,
            rsys: SurfaceRxnSystem = None,
            species: SpeciesManager = None,

            surface: Surface = None,

            time: float = None,
            num_runs: int = None,
            initial_state_seed: int = None,  # TODO(Andrew) should there be separate seeds?
            # the initial state should be different across different runs, right?
            simulation_rng_seed: int = None,
    ):
        # --- Extract input from Spec -----------------------------------------
        rsys = rsys or self.rsys
        species = species or self.sm
        surface = surface or self.surface
        # TODO(Andrew) Allow for initial surface state (see if it's in the Surface class)

        # Collect the names of all the species
        # TODO(Andrew): the SurfaceSpeciesManager should know about the surface
        species_names = rsys.species

        time = time or self.time
        num_runs = num_runs or self.runs  # TODO(Andrew) change name

        initial_state_seed = initial_state_seed or self.initial_state_seed
        simulation_rng_seed = simulation_rng_seed or self.simulation_rng_seed

        # --- Run num_runs simulations ------------------------------------
        runs_init_surfaces = []
        runs_times = []
        runs_species_counts = []
        event_histories = []

        # TODO(Andrew): Some bodges
        rsys.species_manager = self.sm
        self.sm.large_species_dict = {}

        # TODO(Andrew): Replace with a simple dict which I supply
        def site_species_map(sm, specie_name):
            site_species_map = {}
            for species in sm.sm:
                site_species_map[species.name] = species.site

            if specie_name in site_species_map:
                return site_species_map[specie_name]
            else:
                return specie_name
        self.sm.site_species_map = site_species_map

        for run_number in range(num_runs):
            run_rng_seed = simulation_rng_seed + run_number  # TODO(Andrew) this might be bad if someone runs several scrns

            # TODO(Andrew): Allow them to supply an init surface too
            # This makes a dictionary of the structure
            #  SiteName: [list of (species, coverage) pairs]
            coverage_info = collections.defaultdict(list)
            for cov in rsys.by_subclass()[Coverage]:
                coverage_info[self.sm[cov.sm].site].append(
                    (cov.sm, cov.coverage))

            init_surface_state = surface.make_state(
                coverage_info=coverage_info,
                size=self.size,
                # TODO(Andrew): Add seeded RNG
            )

            # TODO(Andrew): bodge
            if surface.structure == 'rectangle':
                init_surface = SquareGrid(*self.size)
                init_surface.set_global_state(init_surface_state)
            elif surface.structure == 'hexagon':
                init_surface = HexGridPlusIntersections(*self.size)
                init_surface.set_global_state(*init_surface_state)
            elif surface.structure == 'voronoi':
                init_surface = copy.copy(self.surface._cg)
                init_surface.set_global_state(init_surface_state)
            else:
                assert False, 'We need an initial surface'

            # TODO(Andrew) it also might not even be used (???)
            # TODO(Andrew) use initial_state_seed here; this should be random (?)
            runs_init_surfaces.append(copy.deepcopy(init_surface))  # TODO(Andrew) best way to track this?

            # Convert input into SurfaceCRNOptionParser format for simulator (sclamons)
            # TODO(Andrew) this is temporary, I should dissect it
            #  as of current, if something is available outside of the manifest, I try to
            #  get it from there. This should ultimately be a method of SCRNSpec
            opts = _get_opts_via_manifest(
                rsys=rsys,
                time=time,
                rng_seed=run_rng_seed
            )

            times, species_counts, event_history = simulate_scrn_single_nodisplay(
                rsys=rsys,  # TODO(Andrew) unpack this so sim can be spec-blind
                species_names=species_names,
                init_surface=init_surface,
                time_max=time,
                rng_seed=run_rng_seed,
                opts=opts,
            )

            runs_times.append(times)
            runs_species_counts.append(species_counts)
            event_histories.append(event_history)

        scts = SurfaceCRNTimeSeries.from_times_and_counts(runs_times, runs_species_counts, self)
        # TODO(Andrew) make this part of the actual SurfaceCRNTimeSeries.__init__:
        scts._event_histories = event_histories
        scts._init_surfaces = runs_init_surfaces
        return scts

    def compare_surface_to_bulk(self):
        scts = self.simulate()

        initial_concs = scts.at(0)
        concs = [Conc(symbol, initial_concs[symbol]) for symbol in initial_concs.index]

        crn = CRNSpec(
            *self.rsys.to_bulk().elements,
            *concs,
            species=self.sm,
            time=self.time,
        )

        cts = crn.simulate()

        return scts, cts

    def to_bulk(self):  # For SurfaceCRNSpec
        rsys = self.rsys.to_bulk()

        # I don't know how to get the initial concentration without dissecting everything first
        #  so I'm putting this on hold
        raise NotImplementedError()