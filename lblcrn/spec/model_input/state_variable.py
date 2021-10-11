from typing import Any, List, DefaultDict, Dict, Optional
from lblcrn.spec.spec_abc import Spec, SpecCollection


class StateVariable(Spec):
    """A class used to store and manage input elements.

    Future work: Study the inheritance relationship from InputElement
    """
    def __init__(self,
                 value: Optional[Any] = None,
                 name: str = "",
                 description: str = "",
                 unit: str = ""
                 ):
        super().__init__(name=name, description=description, value=value)
        self.unit = unit

    def __str__(self):
        """
        Return a human-readable version of the state variable.
        """
        base_str = f"{self.name}={self.value} {self.unit}"
        if self.description:
            return f"{base_str}\nDescription: {self.description}"
        return base_str


class T(StateVariable):
    """A class used to store and manage Temperature.
    """
    def __init__(self,
                 value: Optional[Any] = None,
                 description: str = ""
                 ):
        self.name = "Temperature"
        self.unit = "K"
        super().__init__(value=value, name=self.name, description=description, unit=self.unit)


class P(StateVariable):
    """A class used to store and manage Pressure.
    """
    def __init__(self,
                 value: Optional[Any] = None,
                 description: str = ""
                 ):
        self.name = "Pressure"
        self.unit = "Torr"
        super().__init__(value=value, name=self.name, description=description, unit=self.unit)

    def to_pascal(self) -> float:
        """
        Return the corresponding value in Pascal.
        :return: convert to Pascal.
        """
        # Source: https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-
        # factors/nist-guide-si-appendix-b8, accessed on July 9, 2021. not exact.
        conversion_factor = 1.333224e2
        return conversion_factor * self.value


class GibbsFreeEnergy(StateVariable):
    """A class used to store and manage Pressure.
    """
    def __init__(self,
                 value: Optional[Any] = None,
                 description: str = "",
                 ):
        self.name = "deltaG"
        self.unit = "joules"
        super().__init__(value=value, name=self.name, description=description, unit=self.unit)


# # TODO
# class Surface()
#
#     def __init__(self,
#                  type='fcc111'):
#         """
#         'bcc111'
#         :param type:
#         """
