"""TODO"""

from __future__ import annotations
from typing import Dict, List, Mapping, Set, Tuple, Optional

import copy

import numpy as np
import pandas as pd
import sympy as sym

from dtcs.sim.bulk_crn import bulk_crn
from dtcs.spec.crn.crn_abc import CRNSpecABC
from dtcs.spec.crn.bulk.rxn_system import BulkRxnSystem
from dtcs.spec.species import SpeciesManager
from dtcs.twin.crn import CRNTimeSeries


class BulkCRNSpec(CRNSpecABC):
    """TODO
    """

    _rxn_sys_cls = BulkRxnSystem

    def __init__(self,
                 *elements,
                 time: int = 10,
                 max_step: float = 0.1,
                 rsys: BulkRxnSystem = None,
                 species: SpeciesManager = None,
                 **kwargs):
        super().__init__(*elements,
                         rsys=rsys,
                         species=species,)
        self.time = time
        self.max_step = max_step

    def simulate(
            self: BulkCRNSpec,
            time: Optional[float] = None,
            pressure: Optional[float] = None,
            temperature: Optional[float] = None,
            rates: Optional[Dict] = None,
            end_when_settled: bool = False,
            **options,
        ):
        """Simulate the given reaction system over time.

        Args:
            time: The time until which to simulate.
            options: Forwarded to scipy.integrate.solve_ivp
            TODO

        Returns:
            A CRNTimeSeries object with time series data.
        """
        self.rsys: BulkRxnSystem

        time = time or self.time
        # pressure = pressure or self.pressure
        # temperature = temperature or self.temperature
        rates = rates or self.rsys.get_rates(pressure, temperature)
        # TODO: Throw an error if they ask for rates and temp/pressure
        end_when_settled = end_when_settled or (time is None)

        sol_t, sol_y = bulk_crn.solve_rsys_ode(
            evolution_func=self.rsys.get_evolution_function(rates),
            concentration_func=self.rsys.get_concentration_function(),
            schedule=self.rsys.schedule,
            species=self.rsys.species,
            time_max=time,
            end_when_settled=end_when_settled,
            # max_step=self.max_step,
            **options
        )
        cts = CRNTimeSeries(sol_t, sol_y, copy.deepcopy(self))

        return cts

    def simulate_xps(
            self,
            title: str = "",
            species: List[sym.Symbol] = None,
            ignore: List[sym.Symbol] = None,
            x_range: Optional[np.ndarray] = None,
            scale_factor: float = None,
            experimental: pd.Series = None,
            gas_interval: Tuple[float, float] = None,
            contam_spectra: Optional[Dict[sym.Symbol, pd.Series]] = None,
            deconv_species: Optional[List[sym.Symbol]] = None,
            autoresample: bool = True,
            autoscale: bool = True,
    ):
        cts = self.simulate()
        xps = cts.xps_with(
            t=-1,
            title='',
            species=species,
            ignore=ignore,
            x_range=x_range,
            scale_factor=scale_factor,
            experimental=experimental,
            gas_interval=gas_interval,
            contam_spectra=contam_spectra,
            deconv_species=deconv_species,
            autoresample=autoresample,
            autoscale=autoscale
        )

        return xps
