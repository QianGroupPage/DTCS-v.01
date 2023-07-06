"""This file loads configuration information.

TODO: Make it depend on a yaml file.
"""

from sympy.physics import units

# Whether or not dtcs should warn about depreciated methods
#  'default' will warn, 'ignore' will ignore.
depr_warnings = 'default'

# Default unit system
units = {
    'energy': units.eV,
    'pressure': units.torr,
    'temperature': units.K,
    'time': units.second,
}
ref_tp = {
    'pressure': 0.1,
    'temperature': 273.15,
}
unit_system = list(units.values())
