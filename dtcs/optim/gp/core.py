"""TODO"""
from typing import Iterable

import abc
import logging

_logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
import scipy

try:
    from gpcam.autonomous_experimenter import AutonomousExperimenterGP
except ModuleNotFoundError:
    _logger.info('Didn\'t load module gpcam')

from dtcs.spec.crn.rxn_abc import RevRxnABC
from dtcs.common.util import feature


# --- Classes -----------------------------------------------------------------
class GetRateInstrument:

    def __init__(self, crn, fix_first=False, separate_revrxns=False):
        """Initialize the virtual CRN instrument.

        Args:
            crn: The chemical reaction network.
            fix_first: If True, the first rate will be fixed at the initial
                value. Requires that you call get_init_position at some point.
            separate_revrxns: True, False, or a list of True/False, describing
                which reversible reactions should be considered to have two
                independent rates, or just one rate.
        """
        self.crn = crn

        # Record the first rate to keep it fixed.
        self._fix_first = fix_first
        self._fixed_rate = None

        # Track which reversible reactions should be considered as two separate
        #  rates to gpCAM (alternative to k, 1/k)
        if separate_revrxns is True:
            self._separate_revrxns = [isinstance(rxn, RevRxnABC)
                                      for rxn in crn.rsys.reactions]
        elif separate_revrxns is False:
            self._separate_revrxns = [False] * len(crn.rsys.reactions)
        elif isinstance(separate_revrxns, Iterable):
            if any([separate and not isinstance(rxn, RevRxnABC)
                    for rxn, separate in
                    zip(crn.rsys.reactions, separate_revrxns)]):
                raise ValueError('You specified to separate the rates of a '
                                 'non-reversible reaction.')
            self._separate_revrxns = separate_revrxns
        else:
            raise TypeError(f'Invalid value for separate_revrxns, '
                            f'{separate_revrxns}, please supply a boolean or '
                            f'list of booleans')

        # Debug information
        self._iter_count = 1

        # Score tracking
        self.best_score = float('-inf')
        self.best_rates = None

    def get_init_position(self):
        """Get the initial position (rates) for gpCAM, in a flattened list."""
        rates = self.crn.get_rates()
        flat_rates = []
        for rate, separate in zip(rates, self._separate_revrxns):
            if isinstance(rate, Iterable):
                if separate:
                    flat_rates.extend(rate)
                else:
                    flat_rates.append(rate[0])
            else:
                flat_rates.append(rate)

        # Set the fixed rate
        if self._fix_first:
            self._fixed_rate = flat_rates.pop(0)

        return flat_rates

    def inst_func(self, to_score_list):
        """Evaluate possible rates, using our scoring method.

        Args:
            to_score_list: A list of dicts. The 'position' key is the rates
                that gpCAM wants us to test, and we set the 'value' key with
                our scoring.

        Returns:
            A dictionary matching the structure of rates_to_score, with the
            'value' keys filled in.
        """
        _logger.info(f'Iteration #{self._iter_count}')
        self._iter_count += 1

        for to_score in to_score_list:
            rates = self._pack(to_score['position'])
            _logger.info(f'Rates: {rates}')

            # Score the supplied rates
            score = self.score(rates)
            to_score['value'] = score
            _logger.info(f'Score: {score}')

            # Track if those rates are the best
            if score > self.best_score:
                self.best_score = score
                self.best_rates = rates
                _logger.info(f'Found new best!')

        return to_score_list

    @abc.abstractmethod
    def score(self, rates) -> float:
        """Assign the given rates a score, where positive is better."""
        raise NotImplementedError()

    def _pack(self, position: np.ndarray):
        rates = list(position)
        packed_rates = []

        # Add the first fixed rate manually
        if self._fix_first:
            rates.insert(0, self._fixed_rate)

        for rxn in self._separate_revrxns:
            if rxn:
                packed_rates.append((rates.pop(0), rates.pop(0)))
            else:
                packed_rates.append(rates.pop(0))
        return packed_rates


class DefaultInstrument(GetRateInstrument):
    def __init__(self, crn, experimental, ignore,
                 fix_first=False, separate_revrxns=False):
        super().__init__(crn,
                         fix_first=fix_first,
                         separate_revrxns=separate_revrxns)
        self.exp_env = experimental.to_numpy()
        self.exp_be = experimental.index.to_numpy().astype(float)
        self.ignore = ignore

    def score(self, rates):
        xps, r_squared, rmse = simulate_and_compare(
            crn=self.crn,
            scaled=rates,
            exp_env=self.exp_env,
            exp_be=self.exp_be,
            ignore=self.ignore
        )
        return r_squared


def kernel_l2_single_task(x1, x2, hyperparameters, obj):
    ################################################################
    ###standard anisotropic kernel in an input space with l2########
    ###########################for single task######################
    '''
    x1: 2d numpy array of points
    x2: 2d numpy array of points
    obj: object containing kernel definition

    Return:
    -------
    Kernel Matrix
    '''
    hps = hyperparameters
    distance_matrix = np.zeros((len(x1), len(x2)))
    # you can always get some help:
    # help(obj.squared_exponential_kernel)
    for i in range(len(x1[0]) - 1):
        distance_matrix += abs(np.subtract.outer(x1[:, i], x2[:, i]) / hps[1 + i]) ** 2
    distance_matrix = np.sqrt(distance_matrix)
    # return   hps[0] * obj.squared_exponential_kernel(distance_matrix,1)
    # return   hps[0] * obj.exponential_kernel(distance_matrix,1) + obj.periodic_kernel(abs(x1[:,-1]-x2[:,-1]), hps[-1],hps[-2])
    return hps[0] * obj.matern_kernel_diff1(distance_matrix, 1)
    # return   hps[0] * obj.matern_kernel_diff2(distance_matrix,1)


def simulate_and_compare(crn, scaled, exp_env, exp_be, ignore):
    crn = crn.subs_rates(scaled)
    xps = crn.simulate_xps(title='',
                           experimental=pd.Series(data=exp_env, index=exp_be),
                           ignore=ignore)

    sim_env = np.array([])
    reconv_env = np.array([])

    raw_reconv_env_base = xps.df.deconvoluted.envelope.to_numpy()
    raw_reconv_be_base = xps.df.deconvoluted.index.to_numpy().astype(float)

    raw_sim_env_base = xps.df.simulated.envelope.to_numpy()
    raw_sim_be_base = xps.df.simulated.index.to_numpy().astype(float)

    raw_sim_env, raw_sim_be = np.array([]), np.array([])
    raw_reconv_env, raw_reconv_be = np.array([]), np.array([])

    for i, be in enumerate(raw_sim_be_base):
        if be <= 540:
            raw_sim_env = np.append(raw_sim_env, raw_sim_env_base[i])
            raw_reconv_env = np.append(raw_reconv_env, raw_reconv_env_base[i])
            raw_sim_be = np.append(raw_sim_be, be)
            raw_reconv_be = np.append(raw_reconv_be, be)

    input_i = 0

    # Resample the stored data to match the number of datapoints in the experimental input.
    for j, pt in enumerate(raw_sim_be):
        if pt > exp_be[input_i]:
            if j == 0:
                sim_env = np.append(sim_env, raw_sim_env[j])
                reconv_env = np.append(reconv_env, raw_reconv_env[j])
            else:
                sim_env = np.append(sim_env, raw_sim_env[j - 1])
                reconv_env = np.append(reconv_env, raw_reconv_env[j - 1])
            input_i += 1
        if input_i == len(exp_env):
            break

    if len(sim_env) < len(exp_env):
        sim_env = np.append(sim_env, np.array([sim_env[-1] for _ in range(len(exp_env) - len(sim_env))]))
        reconv_env = np.append(reconv_env, np.array([reconv_env[-1] for _ in range(len(exp_env) - len(reconv_env))]))

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(reconv_env, sim_env)

    rmse = ((exp_env - sim_env) ** 2).mean() ** .5
    r_squared = r_value ** 2

    return xps, r_squared, rmse


@feature('matproj')
def evaluate(crn, inst_cls=DefaultInstrument, inst_args=None, iterations=100):
    instrument = inst_cls(
        crn=crn,
        **inst_args
    )
    init_position = instrument.get_init_position()

    hpbounds = [[0.001, 1000]] + [[0.1, 5000] for _ in init_position]
    exp = AutonomousExperimenterGP(
        parameter_bounds=np.array([(val / 100, val * 100)
                                   for val in init_position]),
        instrument_func=instrument.inst_func,
        hyperparameters=np.ones(len(init_position) + 1),
        hyperparameter_bounds=hpbounds,
        acq_func='ucb',
        x=np.asarray([init_position]),
        #init_dataset_size=num_rates,
        kernel_func=kernel_l2_single_task,
    )

    exp.train()
    exp.go(
        N=iterations,
        acq_func_opt_setting=lambda number: 'global' if number % 2 == 0 else 'global',
   )

    return instrument