from typing import List, Dict, Optional, Tuple
import pandas as pd
import numpy as np
import sympy as sym
from lblcrn import bulk_crn
from lblcrn.experiments import time_series
from lblcrn.crn_sym.rxn_system import RxnSystem


def simulate_xps(rsys: RxnSystem,
                 time: Optional[float] = None,
                 end_when_settled: bool = False,
                 title: str = '',
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
                 **options):
    """Simulate the given reaction system over time.

    Args:
        TODO: add rest
        rsys: ReactionsSystem, the reaction system to simulate
        time: The time until which to simulate.
        species: The species to include in the XPS.
        ignore: The species to not include in the XPS.
        autoresample: Decides if the XPS resamples on edits.
        autoscale: Decides if the XPS will automatically scale
            the gaussians and envelope to match the experimental data.
        experimental: The experimental value of the XPS.
        gas_interval: The interval in which the peak of the gas phase is
            in the XPS.
        scale_factor: The scale factor by which to scale the simulated
            gaussians in the XPS.
        title: The name to give the XPS, used in plotting.
        **options: Forwarded to scipy.integrate.solve_ivp

    Returns:
        An XPSExperiment with the simulation results.
    """
    end_when_settled = end_when_settled or (time is None)

    sol_t, sol_y = bulk_crn.solve_rsys_ode(rsys, time or 0, end_when_settled, **options)
    cts = time_series.CRNTimeSeries(sol_t, sol_y, rsys)

    return cts.xps_with(title=title,
                        species=species,
                        ignore=ignore,
                        x_range=x_range,
                        scale_factor=scale_factor,
                        experimental=experimental,
                        gas_interval=gas_interval,
                        contam_spectra=contam_spectra,
                        deconv_species=deconv_species,
                        autoresample=autoresample,
                        autoscale=autoscale)

def simulate_xps_with_cts(rsys: RxnSystem,
                 time: Optional[float] = None,
                 end_when_settled: bool = False,
                 title: str = '',
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
                 **options):
    """Simulate the given reaction system over time.

    Args:
        TODO: add rest
        rsys: ReactionsSystem, the reaction system to simulate
        time: The time until which to simulate.
        species: The species to include in the XPS.
        ignore: The species to not include in the XPS.
        autoresample: Decides if the XPS resamples on edits.
        autoscale: Decides if the XPS will automatically scale
            the gaussians and envelope to match the experimental data.
        experimental: The experimental value of the XPS.
        gas_interval: The interval in which the peak of the gas phase is
            in the XPS.
        scale_factor: The scale factor by which to scale the simulated
            gaussians in the XPS.
        title: The name to give the XPS, used in plotting.
        **options: Forwarded to scipy.integrate.solve_ivp

    Returns:
        An XPSExperiment with the simulation results as well as a CRNTimeSeries object with the
        time series data.
    """
    end_when_settled = end_when_settled or (time is None)

    sol_t, sol_y = bulk_crn.solve_rsys_ode(rsys, time, end_when_settled, **options)
    cts = time_series.CRNTimeSeries(sol_t, sol_y, rsys)

    return cts.xps_with(title=title,
                        species=species,
                        ignore=ignore,
                        x_range=x_range,
                        scale_factor=scale_factor,
                        experimental=experimental,
                        gas_interval=gas_interval,
                        contam_spectra=contam_spectra,
                        deconv_species=deconv_species,
                        autoresample=autoresample,
                        autoscale=autoscale), cts
