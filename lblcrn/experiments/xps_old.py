
from __future__ import annotations
import copy
from typing import Any, Dict, List, Optional, Tuple, Union
import warnings

from matplotlib import pyplot as plt
import monty.json
import numpy as np
import pandas as pd
from scipy import integrate
from scipy import optimize
from scipy import stats
from sklearn import metrics
import sympy as sym

from lblcrn import bulk_crn
from lblcrn.experiments import experiment
from lblcrn.experiments import time_series
from lblcrn.crn_sym import reaction
from lblcrn.crn_sym import species
from lblcrn.common import util
from lblcrn import _echo

from lblcrn.experiments.xps import XPSObservable

_PLOT_MARGIN = 5
_PLOT_RESOLUTION = 0.001
_SIGMA = 0.75 * np.sqrt(2) / (np.sqrt(2 * np.log(2)) * 2)

_REQUIRED = 'required'
_OPTIONAL = 'optional'

class XPSExperiment:
    # TODO: new methods ------------------------------------------------------
    # TODO: new methods ------------------------------------------------------
    # TODO: new methods ------------------------------------------------------
    # TODO: new methods ------------------------------------------------------
    # TODO: new methods ------------------------------------------------------
    # TODO: new methods ------------------------------------------------------

    @staticmethod
    def _argparse(**inp_args: Any) -> Dict[str, Any]:
        """TODO

        Returns:
            A dict with keys:

                xps: XPSObservable or XPSExperiment
                cache: bool
                overwrite: bool

                species: list of sym.Symbols, perhaps empty
        """
        out_args = {}

        # Convert (xps, df, sm) into xps of type XPSObservable/XPSExperiment
        xps = df = None
        self = inp_args.pop('self', None)
        if isinstance(self, XPSObservable):
            xps = self
        elif isinstance(self, pd.DataFrame):
            xps = df

        xps = inp_args.pop('xps', xps)
        df = inp_args.pop('df', df)
        df = inp_args.pop('dataframe', df)
        species_manager = inp_args.pop('species_manager', None)
        species_manager = inp_args.pop('sm', species_manager)

        if xps is None and df is None:
            raise TypeError('Must supply either xps or df.')
        elif xps is not None and df is not None:
            raise TypeError('Cannot supply both xps and df.')

        if df is not None:
            if species_manager is None:
                raise TypeError('If you supply df, you must also supply sm.')
            xps = XPSObservable(df=df, species_manager=species_manager)

        out_args['xps'] = xps

        # Include 'cache' and 'overwrite' options.
        out_args['cache'] = inp_args.pop('cache', True)
        out_args['overwrite'] = inp_args.pop('overwrite', True)

        # Convert (species, ignore) into 'species'
        out_args['species'] = experiment._get_species_not_ignored(
            species=inp_args.pop('species', None),
            ignore=inp_args.pop('ignore', None),
            all_species=xps.species
        )

        # Look through the df/cache for arguments
        def arg_default(key, default):
            """If out_args[key] is undefined or None:
                (1) Set it to inp_args[key] if inp_args[key] is not None.
                (2) Set it to default.

            Guarantees:
                (1) out_args[key] will not throw an error.
                (2) key will not be in inp_args.
            """
            if key in out_args and out_args[key] is not None:
                return
            value = inp_args.pop(key, None)
            if value is None:
                value = default
            out_args[key] = value

        if isinstance(xps, XPSExperiment):
            # Auto-flags
            arg_default('autoresample', xps.autoresample)
            arg_default('autoscale', xps.autoscale)
            # Cached Inputs
            arg_default('sim_concs', xps.sim_concs)
            arg_default('exp_input', xps.exp_input)
            arg_default('gas_interval', xps.gas_interval)
            arg_default('contam_spectra', xps.contam_spectra)
            arg_default('deconv_species', xps.deconv_species)
            # Maintained Inputs
            arg_default('scale_factor', xps.scale_factor)
            arg_default('contam_concs', xps.contam_concs)
            arg_default('deconv_concs', xps.deconv_concs)
        if xps.df is not None:
            arg_default('x_range', xps.x_range)
            arg_default('gaussians', xps.gaussians)
            arg_default('simulated', xps.simulated)
            arg_default('sim_envelope', xps.sim_envelope)
            arg_default('sim_gaussians', xps.sim_gaussians)
            arg_default('experimental', xps.experimental)
            arg_default('exp_clean', xps.exp_clean)
            arg_default('exp_raw', xps.exp_raw)
            arg_default('contaminants', xps.contaminants)
            arg_default('deconvoluted', xps.deconvoluted)
            arg_default('deconv_envelope', xps.deconv_envelope)
            arg_default('deconv_gaussians', xps.deconv_gaussians)

        # Infer what we can / provide the last defaults
        arg_default('autoresample', False)
        arg_default('autoscale', False)

        arg_default('sim_concs', {})
        arg_default('exp_input', None)
        arg_default('gas_interval', ())
        arg_default('contam_spectra', {})
        arg_default('deconv_species', {})

        arg_default('scale_factor', 1.0)
        arg_default('contam_concs', {})
        arg_default('deconv_concs', {})

        arg_default('x_range', None)
        arg_default('gaussians', None)
        arg_default('simulated', None)
        arg_default('sim_envelope', None)
        arg_default('sim_gaussians', None)
        arg_default('experimental', None)
        arg_default('exp_clean', None)
        arg_default('exp_raw', None)
        arg_default('contaminants', None)
        arg_default('deconvoluted', None)
        arg_default('deconv_envelope', None)
        arg_default('deconv_gaussians', None)

        # At this point, inp_args should be empty. Hence, if there are any
        #  keys still in it, they are unknown.
        if inp_args:
            bad_key = list(inp_args.keys())[0]
            raise TypeError(f'Got an unexpected keyword argument '
                            f'\'{bad_key}\')')

        return out_args

    @staticmethod
    def _cache_and_overwrite(kwargs, **data) -> Union[XPSExperiment,
                                                      XPSObservable,
                                                      pd.DataFrame]:
        assert 'xps' in kwargs
        assert 'cache' in kwargs
        assert 'overwrite' in kwargs

        xps = kwargs.pop('xps')
        cache = kwargs.pop('cache')
        overwrite = kwargs.pop('overwrite')

        if isinstance(xps, XPSExperiment):
            xps.autoscale = data.get('autoresample', xps.autoresample)
            xps.autoscale = data.get('autoscale', xps.autoscale)
            if cache:
                xps._sim_concs = data.get('sim_concs',
                                          xps._sim_concs)
                xps._exp_input = data.get('exp_input',
                                          xps._exp_input)
                xps._gas_interval = data.get('gas_interval',
                                             xps._gas_interval)
                xps._contam_spectra = data.get('contam_spectra',
                                               xps._contam_spectra)
                xps._deconv_species = data.get('deconv_species',
                                               xps._deconv_species)
                xps._scale_factor = data.get('scale_factor',
                                             xps._scale_factor)
                xps._contam_concs = data.get('contam_concs',
                                             xps._contam_concs)
                xps._deconv_concs = data.get('deconv_concs',
                                             xps._deconv_concs)
        xps.df = data.get('df',  # Resample should use this one.
                          xps.df)
        if overwrite and xps.df is not None:
            # TODO: Perhaps I should make my own get function so that I'm not
            #  a bunch of 'dataframe=itself's?
            data.get('gaussians', xps.gaussians)
            data.get('simulated', xps.simulated)
            data.get('sim_envelope', xps.sim_envelope)
            data.get('sim_gaussians', xps.sim_gaussians)
            data.get('experimental', xps.experimental)
            data.get('exp_clean', xps.exp_clean)
            data.get('exp_raw', xps.exp_raw)
            data.get('contaminants', xps.contaminants)
            data.get('deconvoluted', xps.deconvoluted)
            data.get('deconv_envelope', xps.deconv_envelope)
            data.get('deconv_gaussians', xps.deconv_gaussians)

        return xps

    # --- Resample -----------------------------------------------------------
    def resample(self, **kwargs):
        """TODO"""
        kwargs = XPSExperiment._argparse(self=self, **kwargs)
        XPSExperiment._can_resample(**kwargs, complain=True)

        # autoscale = kwargs['autoscale']
        simulate = kwargs['simulate']
        experimnet = kwargs['experiment']
        calc_gas_phase = kwargs['calc_gas_phase']
        decontaminate = kwargs['decontaminate']
        contaminate = kwargs['contaminate']
        deconvolute = kwargs['deconvolute']

        if simulate is None:
            simulate = XPSExperiment._can_simulate(**kwargs)
        elif simulate:
            XPSExperiment._can_simulate(**kwargs, complain=True)

        XPSExperiment._do_resample(**kwargs)

    @staticmethod
    def _can_resample(complain=False, **kwargs) -> bool:
        pass

    @staticmethod
    def _do_resample():
        pass

    # --- Scale --------------------------------------------------------------
    def _scale():
        raise NotImplementedError()

    def _simulate(**kwargs):
        raise NotImplementedError()

    def _simulate_envelope():
        raise NotImplementedError()

    def _simulate_gaussians():
        raise NotImplementedError()

    def _experiment():
        raise NotImplementedError()

    # --- Experimental-Set ---------------------------------------------------
    def exp_set(self: Optional[Union[XPSExperiment, XPSObservable,
                                     pd.DataFrame]] = None,
                exp_input: Optional[pd.Series] = None,
                **kwargs):
        """TODO"""
        kwargs = XPSExperiment._argparse(self=self,
                                         exp_input=exp_input,
                                         **kwargs)
        exp_input = kwargs.pop('exp_input')

        XPSExperiment._can_exp_set(exp_input=exp_input,
                                   complain=True)

        data = XPSExperiment._do_exp_set(exp_input=exp_input)

        XPSExperiment._cache_and_overwrite(kwargs, **data)

    @staticmethod
    def _can_exp_set(exp_input: pd.Series,
                     complain: bool = False) -> bool:
        """TODO"""
        if exp_input is None:
            if complain:
                raise TypeError('No experimental data supplied.')
            return False
        return True

    @staticmethod
    def _do_exp_set(exp_input: pd.Series) -> Dict[str, Any]:
        """TODO"""
        return dict(
            exp_input=exp_input,
            exp_raw=exp_input,
        )

    def _clean():
        raise NotImplementedError()

    def _calc_gas_phase():
        raise NotImplementedError()

    def _contaminants_from():
        raise NotImplementedError()

    def _contam_concs_from():
        raise NotImplementedError()

    def _deconvolute():
        raise NotImplementedError()

    def _deconv_calc_envelope():
        raise NotImplementedError()

    def deconv_calc_gaussians(
            self: Optional[Union[XPSExperiment, XPSObservable]] = None,
            species_manager: Optional[species.SpeciesManager] = None,
            species: Optional[List[sym.Symbol]] = None,
            x_range: Optional[np.ndarray] = None,
            guesses: Optional[Dict[sym.Symbol, float]] = None,
            scale_factor: Optional[float] = None,
            to_fit: Optional[pd.Series] = None,
            **kwargs) -> Union[XPSExperiment, XPSExperiment, pd.DataFrame]:
        """TODO"""
        kwargs.update(dict(
            self=self,
            species_manager=species_manager,
            deconv_species=species,
            x_range=x_range,
            sim_concs=guesses,
            scale_factor=scale_factor,
            exp_clean=to_fit,
        ))
        kwargs = XPSExperiment._argparse(**kwargs)
        xps = kwargs.pop('xps')
        species_manager = xps.species_manager
        species = kwargs.pop('deconv_species')
        x_range = kwargs.pop('x_range')
        guesses = kwargs.pop('sim_concs')
        scale_factor = kwargs.pop('scale_factor')
        to_fit = kwargs.pop('exp_clean')

        XPSExperiment._validate_deconv_calc_gaussians(
            species_manager=species_manager,
            species=species,
            x_range=x_range,
            guesses=guesses,
            scale_factor=scale_factor,
            to_fit=to_fit,
            complain=True
        )

        deconv_gaussians, deconv_concs = XPSExperiment._deconv_calc_gaussians(
            species_manager=species_manager,
            species=species,
            x_range=x_range,
            guesses=guesses,
            scale_factor=scale_factor,
            to_fit=to_fit,
        )

        return XPSExperiment._cache_and_overwrite(
            xps, kwargs,
            deconv_gaussians=deconv_gaussians,
            deconv_species=species,
            deconv_concs=deconv_concs,
        )

    @staticmethod
    def _can_deconv_calc_gaussians(
            species_manager: species.SpeciesManager,
            species: List[sym.Symbol],
            x_range: np.ndarray,
            guesses: Optional[dict],
            scale_factor: float,
            to_fit: float,
            complain: bool = False
    ) -> bool:
        if species_manager is None:
            if complain:
                raise TypeError('Species Manager required to plot gaussians.')
            return False
        if not species:
            if complain:
                raise ValueError('No species supplied for deconvolution.')
            return False
        if x_range is None:
            if complain:
                raise TypeError('X-Range required for deconvolution.')
            return False
        if scale_factor is None:
            if complain:
                raise TypeError('Scale factor required for deconvolution.')
        if to_fit is None:
            if complain:
                raise TypeError('Deconvolution requires data to fit to.')
        return True

    @staticmethod
    def _deconv_calc_gaussians(
            species_manager: species.SpeciesManager,
            species: List[sym.Symbol],
            x_range: np.ndarray,
            guesses: Optional[dict],
            scale_factor: float,
            to_fit: float, ):
        """Get the deconvolution of experimental into species.

        Args:
            species_manager: The species manager to use to give the binding
                energies on which to plot the gaussians.
            species: A list of sym.Symbols, the species for which we're
                guessing concentrations.
            x_range: np.ndarray, the x-range on which to calculate gaussians.
            guesses: A dict {species: conc} of initial guesses.
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
            if specie in guesses:
                conc_guesses[index] = guesses[specie]
            # In case we'd get an out of bounds error.
            conc_guesses[index] = max(0.0, conc_guesses[index])

        # Make the function to pass scipy.optimize.curve_fit. Needs signature
        # f(x-values, *y-values) -> ndarray of size to_fit.size
        def envelope_from_concs(x_range, *concs: float):
            envelope = np.zeros(x_range.size)
            for index, specie in enumerate(species):
                hill = XPSExperiment._get_gaussian(
                    species_manager=species_manager,
                    x_range=x_range,
                    specie=specie,
                    conc=concs[index])
                envelope += hill * scale_factor
            return envelope

        fitted_concs, covariance = optimize.curve_fit(f=envelope_from_concs,
                                                      xdata=x_range,
                                                      ydata=to_fit,
                                                      p0=conc_guesses,
                                                      bounds=(0, np.inf))
        return dict(zip(species, fitted_concs))

