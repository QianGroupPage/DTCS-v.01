"""
TODO(Andrew)
"""
from __future__ import annotations

from typing import Optional

from dtcs.sim import bulk_crn
from dtcs.twin.crn import CRNTimeSeries
from dtcs.common import util

@util.depreciate
def _simulate_bulk_crn(
        rsys: BulkRxnSystem = None,
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
    raise NotImplementedError()
    end_when_settled = end_when_settled or (time is None)

    sol_t, sol_y = bulk_crn.solve_rsys_ode(
        rsys=rsys,
        time_max=time,
        end_when_settled=end_when_settled,
        max_step=crn.max_step,
        **options)
    cts = CRNTimeSeries(sol_t, sol_y, crn)

    return cts
