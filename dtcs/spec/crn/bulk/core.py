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
                 max_step: float = 0.01,
                 rsys: BulkRxnSystem = None,
                 species: SpeciesManager = None,
                 **kwargs):
        super().__init__(*elements,
                         rsys=rsys,
                         species=species,)
        self.time = time
        self.max_step = 0.1

    def simulate(
            self: BulkCRNSpec,
            #rsys: BulkRxnSystem = None,
            time: Optional[float] = None,
            end_when_settled: bool = False,
            #title: str = "",
            #species: List[sym.Symbol] = None,
            #ignore: List[sym.Symbol] = None,
            #x_range: Optional[np.ndarray] = None,
            #scale_factor: float = None,
            #experimental: pd.Series = None,
            #gas_interval: Tuple[float, float] = None,
            #contam_spectra: Optional[Dict[sym.Symbol, pd.Series]] = None,
            #deconv_species: Optional[List[sym.Symbol]] = None,
            #autoresample: bool = True,
            #autoscale: bool = True,
            **options
        ):
        """Simulate the given reaction system over time.

        Args:
            TODO: add rest
            #rsys: ReactionsSystem, the reaction system to simulate
            time: The time until which to simulate.
            species: The species to include in the XPS.
            ignore: The species to not include in the XPS.
            autoresample: Decides if the XPS resamples on edits.
            autoscale: Decides if the XPS will automatically scale the gaussians and envelope to match the experimental data.
            experimental: The experimental value of the XPS.
            gas_interval: The interval in which the peak of the gas phase is in the XPS.
            scale_factor: The scale factor by which to scale the simulated gaussians in the XPS.
            title: The name to give the XPS, used in plotting.
            options: Forwarded to scipy.integrate.solve_ivp

        Returns:
            An XPSExperiment with the simulation results as well as a CRNTimeSeries object with the time series data.
        """
        time = time or self.time
        end_when_settled = end_when_settled or (time is None)

        sol_t, sol_y = bulk_crn.solve_rsys_ode(
            rsys=self.rsys,
            time_max=time,
            end_when_settled=end_when_settled,
            max_step=self.max_step,
            **options)
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
