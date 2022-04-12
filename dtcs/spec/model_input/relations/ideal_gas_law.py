"""
Written by:
Ye Wang, July 9, 2021.
"""
from typing import Optional


class IdealGasLawRelation:
    """
    The Ideal Gas Law related pressure p, volume v,
    """
    # Molecular Gas Constant, https://physics.nist.gov/cgi-bin/cuu/Value?r, accessed on July 9, 2021
    r = 8.314462618  # unit: J mol-1 K-1, standard and relative uncertainties are exact.

    def __init__(self, p_in_torr: float, n: float, t: float, v: Optional[float]=None):
        """
        Initiate an Ideal Gas Law relation. Only volume can be unspecified.

        :param p: total pressure, unit: Torr;
        :param v: volume, unit: L;
        :param n: molecule count, unit: mole;
        :param t: absolute temperature, unit: K.
        """
        self.p = p_in_torr
        self.v = v if v else n * IdealGasLawRelation.r * t / self.torr_to_pascal(self.p)
        self.n = n
        self.t = t
        self.n

    def calculate_n(self, p: float, t: float) -> float:
        """
        Return a new molecule count given a set of new pressure and temperature values.
        :param p: pressure p;
        :param t: temperature t;
        :return: updated n value.
        """
        return self.torr_to_pascal(p) * self.v / (IdealGasLawRelation.r * t)

    @staticmethod
    def torr_to_pascal(value: float) -> float:
        """
        Return the corresponding value in Pascal.

        :value: pressure value in Torr;
        :return: convert to Pascal.
        """
        # Source: https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-
        # factors/nist-guide-si-appendix-b8, accessed on July 9, 2021. not exact.
        conversion_factor = 1.333224e2
        return conversion_factor * value
