from __future__ import annotations

from typing import Union, List, Tuple, Dict

import numpy as np
import pandas as pd
import sympy as sym
from scipy import optimize, stats

from dtcs import _logger

SpeciesName: TypeAlias = sym.Symbol
OrbitalInfo: TypeAlias = List[Tuple[float, float]]


_PLOT_MARGIN = 5
_PLOT_RESOLUTION = 0.001
_SIGMA = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)
_NOISE_MAX = 0.005 # A proportion of the maximum value that will form the upper bound for added noise.


# XPS DataFrame topology
_SIMULATED = 'simulated'
_EXPERIMENTAL = 'experimental'
_CONTAMINANTS = 'contaminants'
_DECONVOLUTED = 'deconvoluted'
_SIM_ENV = (_SIMULATED, 'envelope')
_EXP_CLEAN = (_EXPERIMENTAL, 'clean')
_EXP_RAW = (_EXPERIMENTAL, 'raw')
_GAS_PHASE = (_EXPERIMENTAL, 'gas_phase')
_DECONV_ENV = (_DECONVOLUTED, 'envelope')


def simulate_xps(
        # Input Data
        x_range: Union[np.ndarray, None],
        species: List[SpeciesName],
        orbitals: Dict[str, OrbitalInfo],
        sim_concs: Dict[SpeciesName, float],
        sim_species: List[SpeciesName],
        scale_factor: float,
        exp_data: pd.Series,
        gas_interval: Tuple[float, float],
        contam_spectra: Dict[SpeciesName, pd.Series],
        contam_species: List[SpeciesName],
        deconv_species: List[SpeciesName],
        # Flags
        simulate: bool,
        autoscale: bool,
        experimental: bool,
        gas_phase: bool,
        decontaminate: bool,
        contaminate: bool,
        deconvolute: bool,
        augment: bool,
) -> Tuple[pd.DataFrame,
           float,
           Union[Dict[SpeciesName, float], None],
           Union[Dict[SpeciesName, float], None]]:
    """TODO"""
    # --- Invariants -----------------------------------------------------
    assert not (simulate and not sim_concs)
    assert not (simulate and not sim_species)
    assert not (experimental and exp_data is None)
    assert not (gas_phase and not experimental)
    assert not (gas_phase and (len(gas_interval) > 2 or
                               gas_interval[0] > gas_interval[1]))
    assert not (decontaminate and not experimental)
    assert not (contaminate and not simulate)
    # TODO(Andrew): Fix these assertions
    # assert not ((decontaminate or contaminate) and not contam_spectra)
    # assert not ((decontaminate or contaminate) and not contam_species)
    assert not (contaminate and decontaminate)
    assert not (deconvolute and not experimental)
    assert not (deconvolute and not deconv_species)
    assert scale_factor != 0

    # --- Variables for Later --------------------------------------------
    envelope: Union[pd.Series, None] = None
    deconv_concs: Union[Dict[SpeciesName, float], None] = None
    contam_concs: Union[Dict[SpeciesName, float], None] = None
    sim_cols: List[Tuple[str, SpeciesName]] = []
    deconv_cols: List[Tuple[str, SpeciesName]] = []
    contam_cols: List[Tuple[str, SpeciesName]] = []

    # --- Make X-Range ---------------------------------------------------
    # Find the x_range based on what species we need to see.
    if x_range is None:
        species_to_plot = []
        if simulate:
            species_to_plot.extend(sim_species)
        # Don't bother with deconv_/contam_species, because if they're
        # defined, then the x_range will be exp_data.index anyway.
        x_range = _get_x_range(
            orbitals=orbitals,
            species=species_to_plot,
        )
        _logger.echo(f'Using automatically-generated x-range '
                   f'[{x_range[0]}, ..., {x_range[-1]}]...')

    # --- Make Empty DataFrame -------------------------------------------
    # I add all of the columns early for sorting reasons;
    # comparisons between sym.Symbol and str are not defined, so you can't
    # easily sort the columns.
    columns = []
    if simulate:
        columns.append(_SIM_ENV)
        sim_cols = [(_SIMULATED, specie) for specie in sim_species]
        columns.extend(sim_cols)
    if experimental:
        columns.append(_EXP_RAW)
        columns.append(_EXP_CLEAN)
    # # Gas phase and decontaminate will clean up the raw experimental.
    # if gas_interval or decontaminate:
    if gas_phase:
        columns.append(_GAS_PHASE)
    if decontaminate or contaminate:
        contam_cols = [(_CONTAMINANTS, specie)
                       for specie in contam_species]
    if deconvolute:
        columns.append(_DECONV_ENV)
        deconv_cols = [(_DECONVOLUTED, specie)
                       for specie in deconv_species]
        columns.extend(deconv_cols)

    df = _get_dataframe(x_range=x_range,
                                      columns=columns)

    # --- Calculations ---------------------------------------------------
    if simulate:
        if _logger.do_echo:
            sim_concs_used = {specie: conc for specie, conc in
                              sim_concs.items() if specie in sim_species}
            _logger.echo(f'Simulating XPS experiment for species with '
                         f'concentrations {sim_concs_used}...')
        for sim_col in sim_cols:
            # sim_col is a tuple (self._SIMULATED, specie)
            sim_specie = sim_col[1]
            df[sim_col] = _get_gaussian(
                orbital=orbitals[str(sim_specie)],
                x_range=x_range,
                conc=sim_concs[sim_specie],
            )
        df[_SIM_ENV] = df[sim_cols].sum(axis=1)
        envelope = df[_SIM_ENV]

    if experimental:
        _logger.echo(f'Processing experimental data...')
        df[_EXP_CLEAN] = exp_data
        if _EXP_CLEAN in df:
            df[_EXP_CLEAN] = exp_data.copy()

    if gas_phase:
        _logger.echo(f'Calculating gas phase...')
        gas_gauss = _get_gas_phase(x_range=x_range,
                                       exp_data=exp_data,
                                       gas_interval=gas_interval)
        df[_GAS_PHASE] = gas_gauss
        _logger.echo(f'Removing gas phase from experimental...')
        df[_EXP_CLEAN] -= df[_GAS_PHASE]

    # Must go after simulate
    if autoscale:
        # Andrew: This line here was used for intelligent autoscaling.
        #  species_info = [self.species_manager[specie] for specie in species]
        scale_factor = _get_autoscale(sim_envelope=envelope,
                                          exp_envelope=exp_data)
        _logger.echo(f'Generating auto-scale factor {scale_factor}...')

    # Must go after autoscale
    if contaminate or decontaminate:
        contam_concs = {}
        for cont_col in contam_cols:
            # cont_col is a tuple (self._CONTAMINANT, specie)
            cont_specie = cont_col[1]
            contamination, contam_conc = _get_contamination(
                orbitals=orbitals,
                x_range=x_range,
                contam_specie=cont_specie,
                species=species,
                scale_factor=scale_factor,
                spectrum=contam_spectra[cont_specie]
            )
            contam_concs[cont_specie] = contam_conc
            df[cont_col] = contamination
            _logger.echo(f'Calculated contaminant concentrations at '
                       f'{contam_concs}...')

    if contaminate:
        _logger.echo('Adding contaminants to simulated data...')
        for cont_col in contam_cols:
            df[_SIM_ENV] += df[cont_col]

    # Scale everything which needs to be scaled
    _logger.echo(f'Scaling data by factor {scale_factor}...')
    for sim_col in sim_cols:
        df[sim_col] *= scale_factor
    for cont_col in contam_cols:
        df[cont_col] *= scale_factor
    if _SIM_ENV in df:
        df[_SIM_ENV] *= scale_factor

    # Must happen after scaling
    if decontaminate:
        _logger.echo('Removing contaminants from experimental data...')
        for cont_col in contam_cols:
            df[_EXP_CLEAN] -= df[cont_col]

    # Must go after scaling, gas_phase, and decontaminate
    if deconvolute:
        # Get a dict {specie: conc}, like sim_concs
        deconv_concs = _get_deconvolution(
            x_range=x_range,
            species=deconv_species,
            species_guesses=sim_concs,
            orbitals=orbitals,
            scale_factor=scale_factor,
            to_fit=df[_EXP_CLEAN]
        )
        for dec_col in deconv_cols:
            # dec_col is a tuple (self._DECONV, specie)
            dec_specie = dec_col[1]
            df[dec_col] = _get_gaussian(
                orbital=orbitals[str(dec_specie)],
                x_range=x_range,
                conc=deconv_concs[dec_specie]
            )
        df[_DECONV_ENV] = df[deconv_cols].sum(axis=1)
        _logger.echo(f'Deconvoluted experimental to get concentrations '
                   f'{deconv_concs}...')

        # More scaling
        for dec_col in deconv_cols:
            df[dec_col] *= scale_factor
        if _DECONV_ENV in df:
            df[_DECONV_ENV] *= scale_factor

    if augment:
        m = max(df[_SIM_ENV])
        df[_SIM_ENV] += (_NOISE_MAX * m) * np.random.normal(0, 1, len(df[_SIM_ENV]))

    return df, scale_factor, deconv_concs, contam_concs


def _get_dataframe(x_range: np.ndarray,
                   columns: List) -> pd.DataFrame:
    """Get a properly-formatted dataframe with the given x_range and
    columns."""
    row_index = pd.Index(x_range, name='eV')
    col_index = pd.MultiIndex.from_tuples(columns)
    return pd.DataFrame(data=0, index=row_index, columns=col_index)


def _get_autoscale(sim_envelope: Union[pd.Series, None],
                   exp_envelope: Union[pd.Series, None]) -> float:
    """Gets the factor by which to automatically scale.

    Currently, gives 1.0 unless there's an experimental and a simulated,
    in which case it makes the two envelopes match at the location of the
    highest peak.

    Args:
        df: The pd.DataFrame with the gaussians to scale.
        sim_envelope: The simulated data you're going to scale.
        exp_envelope: Experimental data, if any, to autoscale to.
        species: The species for which the highest peak will be made to
            match.
    """
    scale_factor = 1.0
    if exp_envelope is not None and sim_envelope is not None:
        assert exp_envelope.index.equals(sim_envelope.index), 'Simulated' \
                'and experimental indices should be the same by now.'

        lsq_sol, residuals, _, _ = np.linalg.lstsq(
            a=np.asarray([sim_envelope]).T,  # coefficient matrix, A
            b=exp_envelope,  # goal vector, b
            rcond=-1,
        )  # This should do argmin_x(|Ax-b|^2), which is 1-dimensional least squares (x is a 1-by-1 matrix)

        scale_factor = lsq_sol[0]
        # TODO(Andrew) Maybe add this as a debugging utility, that would
        #  be cool.
        # dtcs._debug_exports = {
        #     'peaks': peaks,
        #     'exp_peaks': exp_envelope.loc[peaks],
        #     'sim_peaks': sim_envelope.loc[peaks],
        # }
    return scale_factor


def _get_contamination(orbitals: Dict[str, OrbitalInfo],
                       x_range: np.ndarray,
                       contam_specie: SpeciesName,
                       species: List[SpeciesName],
                       scale_factor: float,
                       spectrum: pd.Series, ) -> Tuple[np.ndarray, float]:
    """Get the contamination spectrum of the given specie.

    Args:
        x_range: The range you care about decontaminating.
        contam_specie: The specie you want the spetrum/concentration of.
        species: All of the species you want to be involved in the
            deconvolution of the contaminant spectrum. If there are many
            contaminants, put them all here.
        scale_factor: The scale factor of the simulated data. Keeps
            simulated and deconvoluted concentrations consistent. A
            placeholder for an intensity function.
        spectrum: The contaminant spectrum.

    Returns:
        The spectrum (a np.ndarray) and the concentration of the
        contaminant (a float).
    """

    # Deconvolute the spectra to get the conc of this contaminant
    deconved = _get_deconvolution(
        x_range=np.asarray(spectrum.index),
        species=species,
        species_guesses={},
        scale_factor=scale_factor,
        to_fit=spectrum)
    contam_conc = deconved[contam_specie]

    # Make the gaussian of this contaminant on the x_range we care about.
    contam_gauss = _get_gaussian(
        orbital=orbitals[str(contam_specie)],
        x_range=x_range,
        conc=contam_conc
    )
    # TODO: I would put a ratio here, but I don't understand how it should
    # TODO: work, so I'll leave it out for now.

    return contam_gauss, contam_conc


def _get_deconvolution(x_range: np.ndarray,
                       species: List[SpeciesName],
                       species_guesses: Dict[SpeciesName, float],
                       orbitals: Dict[str, OrbitalInfo],
                       scale_factor: float,
                       to_fit: pd.Series) -> Dict[SpeciesName, float]:
    """Get the deconvolution of experimental into species.

    Args:
        x_range: np.ndarray, the x-range on which to calculate gaussians.
        species: A list of species names, the species for which we're
            guessing concentrations.
        species_guesses: A dict {species: conc} of initial guesses.
            This should be sim_concs, as far as I can tell. It need not
            be populated.
        scale_factor: float, the factor by which to scale the gaussians.
            This is so that the concentrations of the simulated and
            deconvoluted data are comparable.
        to_fit: pd.Series, the spectrum to deconvolute.

    Returns:
        A dict of {species: conc} of the guessed concentration.
    """
    # Guess 1 for everything, unless there's simulated data, then use that.
    conc_guesses = np.ones(len(species))
    for index, specie in enumerate(species):
        if specie in species_guesses:
            conc_guesses[index] = species_guesses[specie]
        # In case we'd get an out of bounds error.
        conc_guesses[index] = max(0.0, conc_guesses[index])

    # Make the function to pass scipy.optimize.curve_fit. Needs signature
    # f(x-values, *y-values) -> ndarray of size to_fit.size
    def envelope_from_concs(x_range, *concs: float):
        envelope = np.zeros(x_range.size)
        for index, specie in enumerate(species):
            hill = _get_gaussian(orbital=orbitals[str(specie)],
                                     x_range=x_range,
                                     conc=concs[index])
            envelope += hill * scale_factor
        return envelope

    fitted_concs, covariance = optimize.curve_fit(f=envelope_from_concs,
                                                  xdata=x_range,
                                                  ydata=to_fit,
                                                  p0=conc_guesses,
                                                  bounds=(0, np.inf))
    return dict(zip(species, fitted_concs))


def _get_gas_phase(xps, x_range: np.ndarray,
                   exp_data: pd.Series,
                   gas_interval: Tuple[float, float]) -> np.ndarray:
    """Calculates the gas phase given an experimental.

    It picks as the peak of the gas phase the maximum of exp_data in
    gas_interval.

    Args:
        x_range: np.ndarray, the x-values on which to caluclate the
            gaussian.
        exp_data: pd.Series, the experimental data to use.
        gas_interval: An interval (2-ple), wherein you expect the peak of
            the gas phase to lie. Make it small, as otherwise it might pick
            the foot of a taller gaussian as the 'peak'.

    Returns:
        An np.ndarray of the gas phase gaussian.
    """
    # Get the range to search for the gas peak in.
    search = exp_data[gas_interval[0]:gas_interval[1]]

    # Get the location of the highest part of the experimental data in range
    peak_x = max(search.index, key=lambda index: search[index])
    peak = exp_data[peak_x]

    # Make a gaussian the same height as the experimental gas phase peak
    gas_gaussian = stats.norm.pdf(x_range, peak_x, _SIGMA)
    gas_gaussian *= (peak / max(gas_gaussian))

    return gas_gaussian


def _get_x_range(orbitals: Dict[str, OrbitalInfo],
                 species: List[SpeciesName]) -> np.ndarray:
    """Picks an x-range on which to calculate gaussians.

    Args:
        species: A list of species names, the species you want to be visible
            in the x_range.

    Returns:
        An np.ndarray on which we're going to calculate gaussians.
    """
    binding_energies = []
    for specie in species:
        for orbital in orbitals[str(specie)]:
            binding_energies.append(orbital[0])

    x_lower = min(binding_energies or [0]) - _PLOT_MARGIN
    x_upper = max(binding_energies or [0]) + _PLOT_MARGIN
    x_range = np.arange(x_lower, x_upper, _PLOT_RESOLUTION)

    return x_range


def _get_gaussian(orbital: OrbitalInfo,
                  x_range: np.ndarray,
                  conc: float) -> np.ndarray:
    """Calculates the (unscaled) gaussian for the given species.

    Args:
        species_manager: The species manager to use to get the binding
            energy.
        x_range: np.ndarray, the x-values on which to caluclate the
            gaussian.
        specie: str, the specie for which to get the gaussian.
        conc: float, the concentration of the specie.

    Returns:
        A x_range-sized ndarray with a normal distribtuion centered at the
        binding energy of specie, with st.dev _SIGMA.
    """
    gaussian = np.zeros(x_range.size)

    for binding_energy, splitting in orbital:
        hill = stats.norm.pdf(x_range, binding_energy, _SIGMA)
        gaussian += conc * splitting * hill

    return gaussian

