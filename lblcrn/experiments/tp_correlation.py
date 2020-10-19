import re
from typing import Dict, List, Optional

import numpy as np

boltzmann_constant = 8.61733e-5  # eV/K
kcal_mol_ev_conversion_factor = 23.06035
standard_pressure = 760  # 1 atm in Torr

# The definition of the ZPE
zpe_re = r"^\s*The zero point energy \(ZPE\):\s*([0-9.]+)\s*kcal/mol"

# The start of an environment (situation with a certain temperature) definition
env_re = r"^\s*T\s*=\s*([0-9.]+)\s*K"

# TODO: Clean up this regex
# Match the total line in the environment definition
env_total_row_re = (
    r"\s*total\s*(-?[0-9.]+\s+)(-?[0-9.]+\s+)(-?[0-9.]+\s+)(-?[0-9.]+\s+)(-?[0-9.]+\s+)(-?[0-9.]+)"
)
# The number of values defined in the total line of an environment definition
num_defs_env_total_row = 6
env_total_row_s_column = 2
env_total_row_h_column = 3


class XPSTPCorrelator:
    """Correlate temperature and pressure terms to auto-scale reaction constants.
    """

    def __init__(
        self,
        tp_file: str,
        init_temp: float,
        init_pressure: float,
        adsorption: List[int],
        desorption: List[int],
        init_constants: Optional[List[float]] = None,
        quantum_es: Optional[List[float]] = None,
    ):
        """Create a new XPS Temperature-Pressure Correlator.

        Relevant reactions is a list of indices of reactions for which reaction constants should
        be updated.
        """
        if init_constants is None and quantum_es is None:
            raise ValueError("Either init_constants or quantum_es must be defined")

        self.init_constants = init_constants
        self.quantum_es = quantum_es
        self.init_temp = init_temp
        self.init_pressure = init_pressure
        self.adsorption = adsorption
        self.desorption = desorption
        self.tp_file = tp_file
        self.zpe: float = 0
        self.init_corr_g: float = 0

        # H and S map temperature to H/S values
        self.h: Dict[float, float] = {}
        self.s: Dict[float, float] = {}

        self._parse_tp_file()

    def _corr_g(self, temp: float, pressure: float, h: float, s: float) -> float:
        corr_g_t = (self.zpe + h - (s * temp) / 1000) / kcal_mol_ev_conversion_factor
        corr_g_p = boltzmann_constant * temp * np.log(pressure / standard_pressure)
        return corr_g_t + corr_g_p

    def _parse_tp_file(self) -> None:
        """Parse the temperature and pressure specification jaguar file. Calculating corrG.
        """
        zpe: float = 0.0  # kcal/mol

        f = open(self.tp_file)

        it = (lambda fi: (yield from fi))(f)
        line = next(it, None)

        while line is not None:
            # Zero Point Energy
            zpe_match = re.match(zpe_re, line)
            if zpe_match and len(zpe_match.groups()) == 1:
                zpe = float(zpe_match.groups()[0])

            # Environment definitions
            env_match = re.match(env_re, line)
            if env_match and len(env_match.groups()) == 1:
                temp = float(env_match.groups()[0])
                self._parse_env_definition(it, temp)

            line = next(it, None)

        print("ZPE:", str(zpe))
        print("H:", str(self.h))
        print("S:", str(self.s))

        self.zpe = zpe
        self.init_corr_g = self._corr_g(
            self.init_temp, self.init_pressure, self.h[self.init_temp], self.s[self.init_temp]
        )

    def _parse_env_definition(self, it, temp: float):
        """Given a file iterator at the start of an environment definition, parse the definition.

        Use the given temperature to add H and S values.
        """
        # TODO: Increase the robustness of environment parsing
        # Jump 7 lines
        [next(it, None) for _ in range(7)]

        # Total line
        line = next(it, None)
        if line is None:
            raise ValueError("Invalid tp file format: EOF when parsing an env definition")

        match = re.match(env_total_row_re, line)
        if match is None or len(match.groups()) != num_defs_env_total_row:
            raise ValueError(
                "Invalid tp file format: "
                + str(num_defs_env_total_row)
                + " columns of data not found"
            )

        self.h[temp] = float(match.groups()[env_total_row_h_column])
        self.s[temp] = float(match.groups()[env_total_row_s_column])

    def constant(self, index: int, temp: float, pressure: float) -> float:
        """Given reaction index, temperature, pressure, and the type of reaction, return a reaction constant.
        """
        k_new: float = self.init_constants[index]
        kt = boltzmann_constant * temp
        adsorption = True if index in self.adsorption else False
        desorption = True if index in self.desorption else False

        if self.init_constants:
            corr_g = self._corr_g(temp, pressure, self.h[temp], self.s[temp])
            ki = self.init_constants[index]
            if adsorption:
                k_new = ki * np.exp((corr_g - self.init_corr_g) / kt)
            elif desorption:
                k_new = ki * np.exp((-(corr_g - self.init_corr_g)) / kt)
        else:
            e = self.quantum_es[index]
            if adsorption:
                k_new = np.exp(-(e - self.init_corr_g) / kt)
            elif desorption:
                k_new = np.exp((-(e + self.init_corr_g)) / kt)

        return k_new

    def constants(self, temp: float, pressure: float) -> List[float]:
        """Get all reaction constants for a given temperature and pressure.
        """
        return [self.constant(i, temp, pressure) for i in range(len(self.init_constants))]
